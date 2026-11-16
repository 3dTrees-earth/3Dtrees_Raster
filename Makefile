.PHONY: build test test-raster test-all-extra-bytes test-copc-window-reads test-production-contract test-projected-3d-crs \
	test-parallel-production test-denoising test-original-companion-crs-fallback \
	test-artifact-promotion test-voxel-counts test-copc-all-instance-dimensions \
	build-original-order-copc

IMAGE ?= 3dtrees-raster:test
VCS_REF ?= $(shell git rev-parse HEAD 2>/dev/null || printf unknown)
IMAGE_VERSION ?= local

build:
	docker build \
		--build-arg VCS_REF="$(VCS_REF)" \
		--build-arg IMAGE_VERSION="$(IMAGE_VERSION)" \
		--tag "$(IMAGE)" .

test:
	RASTER_IMAGE="$(IMAGE)" bash tests/test_container.sh

test-raster:
	RASTER_SKIP_BUILD=1 RASTER_IMAGE="$(IMAGE)" bash tests/test_raster.sh

.PHONY: test-production-contract
test-production-contract:
	docker run --rm --network none --cpus 2 --memory 1g --memory-swap 1g \
		--entrypoint Rscript --volume "$(CURDIR)/src:/opt/raster:ro" \
		--volume "$(CURDIR)/tests:/tests:ro" "$(IMAGE)" \
		/tests/test_production_contract.R
	RASTER_IMAGE="$(IMAGE)" \
		bash tests/test_production_contract.sh

.PHONY: test-projected-3d-crs
test-projected-3d-crs:
	RASTER_IMAGE="$(IMAGE)" \
		bash tests/test_projected_3d_crs.sh

.PHONY: test-parallel-production
test-parallel-production:
	RASTER_IMAGE="$(IMAGE)" bash tests/test_parallel_production.sh

.PHONY: test-denoising
test-denoising:
	docker run --rm --network none --cpus 2 --memory 1g --memory-swap 1g \
		--env PYTHONDONTWRITEBYTECODE=1 --entrypoint python3 \
		--volume "$(CURDIR):/workspace:ro" "$(IMAGE)" \
		/workspace/tests/test_denoising_validation.py
	docker run --rm --network none --cpus 2 --memory 4g --memory-swap 4g \
		--entrypoint Rscript --volume "$(CURDIR):/workspace:ro" \
		--workdir /workspace/src "$(IMAGE)" \
		-e 'Sys.setenv(RASTER_SOURCE_ONLY = "1"); source("analysis.R"); source("denoising.R"); source("../tests/test_denoising.R"); source("../tests/test_height_policy.R"); source("../tests/test_segment_catalog_boundaries.R")'
	RASTER_IMAGE="$(IMAGE)" bash tests/test_denoising.sh

test-original-companion-crs-fallback:
	docker run --rm --network none \
		--entrypoint Rscript \
		--volume "$(CURDIR):/repo:ro" \
		"$(IMAGE)" \
		/repo/tests/test_original_companion_crs_fallback.R

test-artifact-promotion:
	docker run --rm --network none \
		--entrypoint Rscript \
		--volume "$(CURDIR):/repo:ro" \
		--workdir /repo/src \
		"$(IMAGE)" \
		-e 'Sys.setenv(RASTER_SOURCE_ONLY = "1"); source("/repo/src/raster.R"); source("/repo/tests/test_artifact_promotion.R"); source("/repo/tests/test_companion_bounds.R"); source("/repo/tests/test_resource_validation.R")'

.PHONY: test-voxel-counts
test-voxel-counts:
	docker run --rm --network none --cpus 2 --memory 4g --memory-swap 4g \
		--env PYTHONDONTWRITEBYTECODE=1 --entrypoint python3 \
		--volume "$(CURDIR):/repo:ro" "$(IMAGE)" \
		/repo/tests/test_voxel_counts.py
	RASTER_IMAGE="$(IMAGE)" \
		bash tests/test_voxel_counts.sh

# Differential fixture acceptance: ordered LAZ and canonical COPC must match for
# PredInstance, PredInstance_SAT, and PredInstance_FM under identical resources.
test-copc-all-instance-dimensions:
	docker run --rm --network none --cpus 2 --memory 1g --memory-swap 1g \
		--env PYTHONDONTWRITEBYTECODE=1 --entrypoint python3 \
		--volume "$(CURDIR):/repo:ro" "$(IMAGE)" \
		/repo/tests/test_indexed_copc_bounds.py
	RASTER_IMAGE="$(IMAGE)" \
	bash tests/test_copc_all_instance_dimensions.sh

# Canonical conversion: add OriginalPointIndex, build spatial COPC with Untwine,
# validate point identity/index completeness, and publish only after validation.
build-original-order-copc:
	RASTER_IMAGE="$(IMAGE)" \
	bash scripts/build_original_order_copc.sh \
		"$(ORIGINAL_LAZ)" \
		"$(COPC_LAZ)"

# Patched rlas: selector 0 reads every ExtraByte, also after the ninth.
.PHONY: test-all-extra-bytes
test-all-extra-bytes:
	docker run --rm --network none --cpus 1 --memory 1g --memory-swap 1g \
		--entrypoint Rscript --volume "$(CURDIR)/tests:/tests:ro" "$(IMAGE)" \
		/tests/test_all_extra_bytes.R

# COPC window reads through the hierarchy equal one full read (no laspy seek).
# COPC=/path/to/cloud.copc.laz [WINDOW_M=20] [AREA_M=60]
.PHONY: test-copc-window-reads
test-copc-window-reads:
	docker run --rm --network none --cpus 2 --memory 8g --memory-swap 8g \
		--env PYTHONDONTWRITEBYTECODE=1 --entrypoint python3 \
		--volume "$(CURDIR)/tests:/tests:ro" --volume "$(COPC):/in/cloud.copc.laz:ro" "$(IMAGE)" \
		/tests/compare_copc_window_reads.py /in/cloud.copc.laz $(or $(WINDOW_M),20) $(or $(AREA_M),60)
