.PHONY: build build-julia-memory-safe test test-julia-memory-safe test-original-companion-crs-fallback test-artifact-promotion test-copc-all-instance-dimensions test-dataset150 test-julia-memory-safe-dataset150 test-julia-memory-safe-dataset150-operational-aoi build-original-order-copc test-valid-updated-copc-alignment test-valid-updated-cohort test-valid-updated-oracle-cohort

IMAGE ?= 3dtrees-foreststructure:local
VCS_REF ?= $(shell git rev-parse HEAD 2>/dev/null || printf unknown)
IMAGE_VERSION ?= local
JULIA_MEMORY_SAFE_IMAGE ?= 3dtrees-foreststructure:copc-local
JULIA_ORIGINAL_SCRIPT ?= reference/Indices_Final_run.R

build:
	docker build \
		--build-arg VCS_REF="$(VCS_REF)" \
		--build-arg IMAGE_VERSION="$(IMAGE_VERSION)" \
		--tag "$(IMAGE)" .

build-julia-memory-safe:
	docker build \
		--file Dockerfile.julia-memory-safe \
		--build-arg VCS_REF="$(VCS_REF)" \
		--build-arg IMAGE_VERSION="julia-memory-safe-$(IMAGE_VERSION)" \
		--tag "$(JULIA_MEMORY_SAFE_IMAGE)" .

test:
	bash tests/test_container.sh

.PHONY: test-galaxy-resource-request
test-galaxy-resource-request:
	PYTHONDONTWRITEBYTECODE=1 python3 tests/test_galaxy_resource_request.py

test-julia-memory-safe:
	FORESTSTRUCTURE_JULIA_IMAGE="$(JULIA_MEMORY_SAFE_IMAGE)" \
	bash tests/test_julia_memory_safe.sh

.PHONY: test-production-contract
test-production-contract:
	docker run --rm --network none --cpus 2 --memory 1g --memory-swap 1g \
		--entrypoint Rscript --volume "$(CURDIR)/src:/opt/foreststructure:ro" \
		--volume "$(CURDIR)/tests:/tests:ro" "$(JULIA_MEMORY_SAFE_IMAGE)" \
		/tests/test_production_contract.R
	FORESTSTRUCTURE_JULIA_IMAGE="$(JULIA_MEMORY_SAFE_IMAGE)" \
		bash tests/test_production_contract.sh

.PHONY: test-projected-3d-crs
test-projected-3d-crs:
	FORESTSTRUCTURE_JULIA_IMAGE="$(JULIA_MEMORY_SAFE_IMAGE)" \
		bash tests/test_projected_3d_crs.sh

.PHONY: test-parallel-production
test-parallel-production:
	bash tests/test_parallel_production.sh

.PHONY: test-denoising
test-denoising:
	docker run --rm --network none --cpus 2 --memory 1g --memory-swap 1g \
		--env PYTHONDONTWRITEBYTECODE=1 --entrypoint python3 \
		--volume "$(CURDIR):/workspace:ro" "$(JULIA_MEMORY_SAFE_IMAGE)" \
		/workspace/tests/test_denoising_validation.py
	docker run --rm --network none --cpus 2 --memory 4g --memory-swap 4g \
		--entrypoint Rscript --volume "$(CURDIR):/workspace:ro" \
		--workdir /workspace/src "$(JULIA_MEMORY_SAFE_IMAGE)" \
		-e 'Sys.setenv(FORESTSTRUCTURE_SOURCE_ONLY = "1"); source("run.R"); source("denoising.R"); source("../tests/test_denoising.R"); source("../tests/test_segment_catalog_boundaries.R")'
	bash tests/test_denoising.sh

test-original-companion-crs-fallback:
	docker run --rm --network none \
		--entrypoint Rscript \
		--volume "$(CURDIR):/repo:ro" \
		"$(JULIA_MEMORY_SAFE_IMAGE)" \
		/repo/tests/test_original_companion_crs_fallback.R

test-artifact-promotion:
	docker run --rm --network none \
		--entrypoint Rscript \
		--volume "$(CURDIR):/repo:ro" \
		--workdir /repo/src \
		"$(JULIA_MEMORY_SAFE_IMAGE)" \
		-e 'Sys.setenv(FORESTSTRUCTURE_SOURCE_ONLY = "1"); source("/repo/src/run_julia_memory_safe.R"); source("/repo/tests/test_artifact_promotion.R"); source("/repo/tests/test_companion_bounds.R"); source("/repo/tests/test_resource_validation.R")'

.PHONY: test-voxel-counts
test-voxel-counts:
	docker run --rm --network none --cpus 2 --memory 4g --memory-swap 4g \
		--env PYTHONDONTWRITEBYTECODE=1 --entrypoint python3 \
		--volume "$(CURDIR):/repo:ro" "$(JULIA_MEMORY_SAFE_IMAGE)" \
		/repo/tests/test_voxel_counts.py
	FORESTSTRUCTURE_JULIA_IMAGE="$(JULIA_MEMORY_SAFE_IMAGE)" \
		bash tests/test_voxel_counts.sh

# Differential fixture acceptance: ordered LAZ and canonical COPC must match for
# PredInstance, PredInstance_SAT, and PredInstance_FM under identical resources.
test-copc-all-instance-dimensions:
	docker run --rm --network none --cpus 2 --memory 1g --memory-swap 1g \
		--env PYTHONDONTWRITEBYTECODE=1 --entrypoint python3 \
		--volume "$(CURDIR):/repo:ro" "$(JULIA_MEMORY_SAFE_IMAGE)" \
		/repo/tests/test_indexed_copc_bounds.py
	FORESTSTRUCTURE_JULIA_IMAGE="$(JULIA_MEMORY_SAFE_IMAGE)" \
	bash tests/test_copc_all_instance_dimensions.sh

test-dataset150:
	bash tests/test_dataset150.sh "$(DATASET150_LAZ)"

test-julia-memory-safe-dataset150:
	FORESTSTRUCTURE_JULIA_IMAGE="$(JULIA_MEMORY_SAFE_IMAGE)" \
	bash tests/test_julia_memory_safe_dataset150.sh \
		"$(DATASET150_LAZ)" \
		"$(DATASET150_GPKG)" \
		"$(JULIA_ORIGINAL_SCRIPT)"

test-julia-memory-safe-dataset150-operational-aoi:
	FORESTSTRUCTURE_JULIA_IMAGE="$(JULIA_MEMORY_SAFE_IMAGE)" \
	bash tests/test_julia_memory_safe_dataset150_operational_aoi.sh \
		"$(DATASET150_LAZ)" \
		"$(DATASET150_OPERATIONAL_AOI)" \
		"$(JULIA_ORIGINAL_SCRIPT)"

# Dataset acceptance: run canonical COPC, prove COPC provenance, then compare all
# scientific tables/vectors/rasters against the supplied Julia-faithful oracle.
test-valid-updated-copc-alignment:
	FORESTSTRUCTURE_JULIA_IMAGE="$(JULIA_MEMORY_SAFE_IMAGE)" \
	bash tests/test_valid_updated_copc_alignment.sh \
		"$(DATASET_ID)" \
		"$(COPC_LAZ)" \
		"$(AOI_GEOJSON)" \
		"$(VALID_UPDATED_DIR)" \
		"$(INSTANCE_DIMENSION)" \
		"$(ORIGINAL_LAZ)"

# Canonical conversion: add OriginalPointIndex, build spatial COPC with Untwine,
# validate point identity/index completeness, and publish only after validation.
build-original-order-copc:
	FORESTSTRUCTURE_JULIA_IMAGE="$(JULIA_MEMORY_SAFE_IMAGE)" \
	bash tests/build_original_order_copc.sh \
		"$(ORIGINAL_LAZ)" \
		"$(COPC_LAZ)"

test-valid-updated-cohort:
	FORESTSTRUCTURE_JULIA_IMAGE="$(JULIA_MEMORY_SAFE_IMAGE)" \
	bash tests/test_valid_updated_cohort.sh \
		"$(or $(COHORT_MANIFEST),tests/fixtures/valid_updated_cohort15.tsv)"

test-valid-updated-oracle-cohort:
	FORESTSTRUCTURE_JULIA_IMAGE="$(JULIA_MEMORY_SAFE_IMAGE)" \
	bash tests/test_valid_updated_oracle_cohort.sh \
		"$(or $(COHORT_MANIFEST),tests/fixtures/valid_updated_cohort20_new.tsv)"
