FROM rocker/geospatial:4.5

# Native BLAS pools otherwise start one thread per host core in every worker.
# Parallelism belongs to the bounded catalog/tile process scheduler.
ENV OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1

ARG LIDR_VERSION=4.3.2
ARG LIDR_SHA256=3df920ff4146b0c29680c68e1fd654091906addbbdf9e3ae7e2fa8e964fec4e0

RUN install2.r --error --skipinstalled \
    argparse \
    BH \
    classInt \
    data.table \
    future \
    geometry \
    glue \
    lazyeval \
    parallelly \
    Rcpp \
    RcppArmadillo \
    rgl \
    rlas \
    sf \
    stars \
    terra \
    && wget --quiet \
      "https://cran.r-project.org/src/contrib/Archive/lidR/lidR_${LIDR_VERSION}.tar.gz" \
      -O /tmp/lidR.tar.gz \
    && printf '%s  %s\n' "${LIDR_SHA256}" /tmp/lidR.tar.gz \
      | sha256sum --check --strict \
    && R CMD INSTALL /tmp/lidR.tar.gz \
    && rm /tmp/lidR.tar.gz

ARG VCS_REF=unknown
ARG IMAGE_VERSION=local

# Same streaming reader as the memory-safe runner; GDAL supplies block raster I/O.
RUN apt-get update \
    && apt-get install --yes --no-install-recommends python3-pip python3-gdal \
    && rm -rf /var/lib/apt/lists/* \
    && python3 -m pip install --no-cache-dir --break-system-packages \
      laspy==2.6.1 lazrs==0.7.0 pyproj==3.7.2 \
    && python3 -c 'import laspy, lazrs, pyproj; from osgeo import gdal, ogr, osr'

LABEL org.opencontainers.image.source="https://github.com/3dTrees-earth/3dtrees_foreststructure" \
      org.opencontainers.image.description="Audit-aware forest-structure metrics for LAS/LAZ point clouds" \
      org.opencontainers.image.licenses="GPL-3.0" \
      org.opencontainers.image.revision="${VCS_REF}" \
      org.opencontainers.image.version="${IMAGE_VERSION}" \
      earth.3dtrees.memory-allocation.max-bytes="100000000000" \
      earth.3dtrees.memory-headroom.fraction="0.15" \
      earth.3dtrees.cpu.max="10" \
      earth.3dtrees.algorithm-base.revision="9d1fba9cb14adac51336ee53b4f5f5535395331a" \
      earth.3dtrees.dtm-tin-retry="legacy-floating-point" \
      earth.3dtrees.instance-dimensions="PredInstance,PredInstance_SAT,PredInstance_FM"

COPY src /opt/foreststructure
RUN chmod -R a+rX /opt/foreststructure \
    && chmod a+rx /opt/foreststructure/container_entrypoint.sh

ENTRYPOINT ["/opt/foreststructure/container_entrypoint.sh"]
