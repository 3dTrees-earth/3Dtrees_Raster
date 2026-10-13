"""Check denoising acceptance without numerical tolerances for execution parity."""
from __future__ import annotations

import argparse
import csv
import filecmp
import hashlib
import json
import math
from pathlib import Path

import laspy
import numpy as np
from osgeo import gdal

from compare_parallel_production import compare_raster, read_json, science_files


def csv_rows(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def compare_products(expected: Path, actual: Path) -> int:
    names = {name for name in science_files(expected) if not name.endswith("_noise_points.csv")}
    candidate_names = {name for name in science_files(actual) if not name.endswith("_noise_points.csv")}
    assert names == candidate_names, (expected, actual, names ^ candidate_names)
    for name in sorted(names):
        left, right = expected / name, actual / name
        if left.suffix == ".tif":
            compare_raster(left, right)
        elif left.suffix == ".geojson":
            assert read_json(left) == read_json(right), right
        else:
            assert filecmp.cmp(left, right, shallow=False), right
    return len(names)


def maximum(path: Path) -> float:
    raster = gdal.Open(str(path))
    band = raster.GetRasterBand(1)
    values = band.ReadAsArray()
    valid = band.GetMaskBand().ReadAsArray() != 0
    return float(np.max(values[valid]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    root = parser.parse_args().root
    gdal.UseExceptions()
    gdal.SetConfigOption("GDAL_PAM_ENABLED", "NO")
    runtime = read_json(root / "runtime.json")
    fixture = read_json(root / "input/fixture.json")
    assert fixture["fixture_version"] == 2
    source_hash = hashlib.sha256((root / "input/input.copc.laz").read_bytes()).hexdigest()
    assert source_hash == runtime["source_sha256"], "Input cloud changed during analysis"
    with laspy.open(root / "input/input.copc.laz") as reader:
        assert reader.header.point_count == fixture["points"]
        extras = list(reader.header.point_format.extra_dimension_names)
        assert extras.index("PredInstance") == 14
        assert extras.index("OriginalPointIndex") > 8
    for dimension in fixture["supported_instance_ids"]:
        order = read_json(root / f"input/order_{dimension}.json")
        assert order["multiset_equal"] and order["original_point_index_unique"], order
        assert order["original_point_index_restores_expected_order"], order
    outputs = {name: root / name / "output" for name in ("none", "sor-w1", "sor-w2")}
    resources = {}
    height_domains = {}
    observed_height_maxima = {}
    supported_lower, supported_upper = fixture["expected_supported_height_m"]
    extreme_lower = fixture["expected_high_outlier_minimum_height_m"]
    for name, output in outputs.items():
        provenance = read_json(output / "foreststructure_julia_memory_safe_run.json")
        assert provenance["source_point_cloud_sha256"] == source_hash
        performance = csv_rows(output / "foreststructure_performance.csv")[0]
        scientific = provenance["scientific_parameters"]
        exact_text = performance["observed_maximum_voxel_height_m_exact"]
        assert isinstance(scientific["observed_maximum_voxel_height_m_exact"], str)
        assert scientific["observed_maximum_voxel_height_m_exact"] == exact_text
        observed_maximum = float(exact_text)
        expected_domain = max(1, math.ceil(observed_maximum))
        assert math.isfinite(observed_maximum) and observed_maximum >= 0
        assert math.isfinite(scientific["observed_maximum_voxel_height_m"])
        assert scientific["observed_maximum_voxel_height_m"] <= expected_domain
        assert scientific["maximum_height"] is None
        assert scientific["normalized_height_policy"] == "finite_nonnegative_without_upper_cutoff"
        assert scientific["metric_height_policy"] == performance["metric_height_policy"] == "observed_voxel_height_v1"
        assert scientific["metric_height_upper_m"] == float(performance["metric_height_upper_m"]) == expected_domain
        assert float(performance["maximum_height"]) == math.inf
        for dimension in fixture["supported_instance_ids"]:
            tile_rows = csv_rows(output / f"foreststructure_{dimension}_results.csv")
            assert tile_rows and {float(row["metric_height_upper_m"]) for row in tile_rows} == {expected_domain}
        height_domains[name] = expected_domain
        observed_height_maxima[name] = observed_maximum
        resources[name] = dict(
            parent_peak_rss_mib=provenance["parent_peak_rss_mib"],
            cgroup_memory_peak_mib=provenance["cgroup_memory_peak_mib"],
            temporary_disk_peak_mib=provenance["temporary_disk_peak_mib"],
            analysis_seconds=float(performance["total_seconds"]),
            denoise_seconds=float(performance.get("denoise_seconds", 0)),
            voxel_seconds=read_json(output / "foreststructure_voxel_counts.json")["elapsed_seconds"],
        )
    assert not (outputs["none"] / "foreststructure_noise_points.csv").exists()
    assert not (outputs["none"] / "foreststructure_denoising.json").exists()
    masks = []
    for name in ("sor-w1", "sor-w2"):
        output = outputs[name]
        rows = csv_rows(output / "foreststructure_noise_points.csv")
        assert rows and set(rows[0]) == {"OriginalPointIndex"}
        mask = [int(row["OriginalPointIndex"]) for row in rows]
        assert len(mask) == len(set(mask)), "Noise mask contains duplicate source indices"
        assert all(0 <= value < fixture["points"] for value in mask)
        assert fixture["high_outlier"]["original_point_index"] in mask
        assert fixture["low_outlier"]["original_point_index"] in mask
        assert len(mask) == 2, "Conservative default flagged ordinary synthetic returns"
        assert not set(mask).intersection(fixture["supported_crown_indices"]), "Supported tall crown was removed"
        masks.append(sorted(mask))
        metadata = read_json(output / "foreststructure_denoising.json")
        assert metadata["point_count"] == fixture["points"]
        assert metadata["rejected_points"] == len(mask)
        assert metadata["retained_points"] + len(mask) == fixture["points"]
        assert metadata["raw_voxel_counts_filtered"] is False
        assert metadata["coordinate_changes"] is False
        assert metadata["k"] == 10 and metadata["m"] == 12
        assert metadata["window_size_m"] == 60 and metadata["buffer_m"] == 20
        assert metadata["partition_count"] > 1
        assert metadata["enabled"] is True
        assert metadata["algorithm"] == "lidR::sor"
        assert metadata["method"] == "fixed_window_local_sor"
        assert sum(part["core_points"] for part in metadata["partitions"]) == fixture["points"]
        provenance = read_json(output / "foreststructure_julia_memory_safe_run.json")
        assert provenance["point_order_dimension"] == "OriginalPointIndex"
        assert set(provenance["processed_instance_dimensions"]) == set(fixture["supported_instance_ids"])
        assert provenance["catalog_workers"] == (1 if name == "sor-w1" else 2)
        assert provenance["mode"] == "julia-memory-safe-local-sor"
        assert provenance["denoising"]["source_sha256"] == source_hash
        denoising_plans = [plan for plan in provenance["resource_admission"] if plan["phase"] == "denoise"]
        assert len(denoising_plans) == 1
        admitted_workers = max(wave["workers"] for wave in denoising_plans[0]["waves"])
        assert admitted_workers == provenance["catalog_workers"]
        resources[name]["maximum_admitted_denoise_workers"] = admitted_workers
        voxel = read_json(output / "foreststructure_voxel_counts.json")
        assert voxel["point_accounting"]["read"] == fixture["points"]
        assert voxel["dtm"]["denoising"] == provenance["denoising"]
        assert voxel["noise_mask_applied_to_counts"] is False
        # Raw export retains the high point rejected by the science mask. It is
        # normalized by this run's terrain, so compare its exact resulting bin.
        x, y, z = fixture["high_outlier"]["xyz"]
        terrain = gdal.Open(str(output / "foreststructure_dtm.tif"))
        transform = terrain.GetGeoTransform()
        col, row = int((x - transform[0]) / transform[1]), int((y - transform[3]) / transform[5])
        elevation = float(terrain.GetRasterBand(1).ReadAsArray(col, row, 1, 1)[0, 0])
        height_lower = np.floor((z - elevation) / 0.5) * 0.5
        matches = [record for record in csv_rows(output / "foreststructure_voxel_counts.csv")
                   if float(record["x_center"]) == np.floor(x / 0.5) * 0.5 + 0.25
                   and float(record["y_center"]) == np.floor(y / 0.5) * 0.5 + 0.25
                   and float(record["height_lower_m"]) == height_lower]
        assert len(matches) == 1 and int(matches[0]["point_count"]) == 1, (name, matches)
        assert supported_lower < maximum(output / "foreststructure_chm.tif") < supported_upper
        for dimension, target_id in fixture["supported_instance_ids"].items():
            diagnostics_name = f"foreststructure_{dimension}_segment_diagnostics.csv"
            before = csv_rows(outputs["none"] / diagnostics_name)
            after = csv_rows(output / diagnostics_name)
            # The established diagnostic schema names this column PredInstance
            # for each separately named source-dimension product.
            assert {record["PredInstance"] for record in before} == {record["PredInstance"] for record in after}
            selected_before = next(record for record in before if int(record["PredInstance"]) == target_id)
            selected_after = next(record for record in after if int(record["PredInstance"]) == target_id)
            assert float(selected_before["apex_z"]) > extreme_lower
            assert supported_lower < float(selected_after["apex_z"]) < supported_upper
    assert masks[0] == masks[1], "Worker/adaptive execution changed scientific rejection mask"
    summaries = [read_json(outputs[name] / "foreststructure_denoising.json")
                 for name in ("sor-w1", "sor-w2")]
    for summary in summaries:
        assert summary.pop("elapsed_seconds") >= 0
    assert summaries[0] == summaries[1]
    assert maximum(outputs["none"] / "foreststructure_chm.tif") > extreme_lower
    worker_count = compare_products(outputs["sor-w1"], outputs["sor-w2"])
    report = dict(status="valid", source_sha256=source_hash, input_points=fixture["points"],
        rejected_points=len(masks[0]), supported_crown_points=len(fixture["supported_crown_indices"]),
        injected_outliers_rejected=2, additional_source_points_flagged=len(masks[0]) - 2,
        denoising_parameters={key: summaries[0][key] for key in (
            "algorithm", "method", "lidR_version", "k", "m", "window_size_m", "buffer_m", "origin")},
        default_cli_enables_sor=True,
        metric_height_policy="observed_voxel_height_v1", metric_height_upper_m=height_domains,
        observed_maximum_voxel_height_m=observed_height_maxima,
        supported_crown_height_interval_m=fixture["expected_supported_height_m"],
        execution_invariance_products_compared=worker_count,
        raw_export_retains_rejected_high_point=True,
        elapsed_seconds={name: int((root / name / "elapsed_seconds.txt").read_text()) for name in outputs},
        resources=resources,
        resource_measurement_scope="Runner-reported peaks at end of analysis, before final validation/promotion; not an independently sampled whole-container peak",
        runtime=runtime)
    (root / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Denoising acceptance passed: {root / 'comparison.json'}")


if __name__ == "__main__":
    main()
