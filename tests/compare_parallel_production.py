"""Exact scientific comparison for the small parallel-production fixture.

Performance/provenance records necessarily differ and are checked separately;
no scientific field, CSV row order, raster value or GeoJSON feature is ignored.
Raster storage layout may differ, but CRS, grid, bands and metadata must agree.
"""
from __future__ import annotations

import argparse
import filecmp
import hashlib
import json
from pathlib import Path

import laspy
import numpy as np
from osgeo import gdal


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def compare_raster(expected: Path, actual: Path) -> dict:
    left, right = gdal.Open(str(expected)), gdal.Open(str(actual))
    assert (left.RasterXSize, left.RasterYSize, left.RasterCount) == (
        right.RasterXSize, right.RasterYSize, right.RasterCount
    ), actual
    assert left.GetGeoTransform() == right.GetGeoTransform(), actual
    assert left.GetSpatialRef() and right.GetSpatialRef(), actual
    assert left.GetSpatialRef().IsSame(right.GetSpatialRef()), actual
    assert left.GetMetadata() == right.GetMetadata(), actual
    for index in range(1, left.RasterCount + 1):
        source, candidate = left.GetRasterBand(index), right.GetRasterBand(index)
        assert source.DataType == candidate.DataType, (actual, index)
        assert source.GetDescription() == candidate.GetDescription(), (actual, index)
        assert source.GetMetadata() == candidate.GetMetadata(), (actual, index)
        assert source.GetUnitType() == candidate.GetUnitType(), (actual, index)
        assert source.GetScale() == candidate.GetScale(), (actual, index)
        assert source.GetOffset() == candidate.GetOffset(), (actual, index)
        np.testing.assert_equal(source.GetNoDataValue(), candidate.GetNoDataValue())
        for y in range(0, left.RasterYSize, 128):
            for x in range(0, left.RasterXSize, 128):
                width, height = min(128, left.RasterXSize - x), min(128, left.RasterYSize - y)
                np.testing.assert_array_equal(
                    source.ReadAsArray(x, y, width, height),
                    candidate.ReadAsArray(x, y, width, height),
                    err_msg=f'{actual.name}, band={index}, x={x}, y={y}',
                )
                np.testing.assert_array_equal(
                    source.GetMaskBand().ReadAsArray(x, y, width, height),
                    candidate.GetMaskBand().ReadAsArray(x, y, width, height),
                )
    return dict(status='valid', bands=left.RasterCount,
                width=left.RasterXSize, height=left.RasterYSize)


def science_files(directory: Path) -> set[str]:
    return {path.name for path in directory.iterdir()
            if path.suffix in ('.csv', '.geojson', '.tif') and 'performance' not in path.name}


def compare_case(root: Path, scope: str, workers: int, point_count: int) -> dict:
    expected = root / f'{scope}-baseline-w1' / 'output'
    case = root / f'{scope}-candidate-w{workers}'
    actual = case / 'output'
    names = science_files(expected)
    assert names == science_files(actual), case
    assert len(names) == 14, (case, sorted(names))  # 3 label sets, 4 rasters, voxel CSV.
    checks = {}
    for name in sorted(names):
        source, candidate = expected / name, actual / name
        if source.suffix == '.csv':
            assert filecmp.cmp(source, candidate, shallow=False), candidate
            checks[name] = dict(status='valid', comparison='byte equality')
        elif source.suffix == '.geojson':
            assert read_json(source) == read_json(candidate), candidate
            checks[name] = dict(status='valid', comparison='complete parsed equality')
        else:
            checks[name] = compare_raster(source, candidate)
    baseline = read_json(expected / 'foreststructure_julia_memory_safe_run.json')
    provenance = read_json(actual / 'foreststructure_julia_memory_safe_run.json')
    for field in ('status', 'source_point_cloud_sha256', 'point_order_dimension',
                  'processed_instance_dimensions', 'scientific_parameters',
                  'source_crs_wkt', 'raster_crs_wkt', 'footprint_source'):
        assert baseline[field] == provenance[field], (case, field)
    assert provenance['catalog_workers'] == workers, case
    assert provenance['processed_instance_dimensions'] == [
        'PredInstance', 'PredInstance_SAT', 'PredInstance_FM'
    ], case
    assert provenance['point_order_dimension'] == 'OriginalPointIndex', case
    groups = provenance['instance_read_groups']
    assert len(groups) == 1, (case, groups)
    assert groups[0]['dimensions'] == provenance['processed_instance_dimensions'], (case, groups)
    admission = provenance['resource_admission']
    assert {plan['phase'] for plan in admission} == {'DTM', 'CHM', 'segment', 'tiles'}, case
    admitted_workers = []
    for plan in admission:
        for wave in plan['waves']:
            admitted_workers.append(wave['workers'])
            assert 1 <= wave['workers'] <= workers, (case, wave)
            assert (wave['reserved_bytes'] + wave['memory_before_bytes'] +
                    plan['reserve_bytes'] <= plan['budget_bytes'] + 1), (case, wave)
    assert max(admitted_workers) == workers, (case, admitted_workers)
    assert provenance['footprint_source'] == (
        'point_cloud_extent' if scope == 'full' else baseline['footprint_source']
    ), case
    expected_counts = read_json(expected / 'foreststructure_voxel_counts.json')
    counts = read_json(actual / 'foreststructure_voxel_counts.json')
    # Paths, hashes of raster encodings and resource estimates can legitimately
    # differ. The voxel semantics, accounting and geometry cannot.
    for field in ('metric', 'units', 'voxel_size_m', 'point_accounting',
                  'domain_cells', 'occupied_voxels', 'grid', 'band_count',
                  'height_reference', 'height_interval', 'band_bounds_formula',
                  'normalization', 'negative_height_policy', 'outlier_policy',
                  'boundary_policy', 'zero', 'nodata', 'validity_mask'):
        assert expected_counts[field] == counts[field], (case, field)
    assert counts['point_accounting']['read'] == point_count, case
    if scope == 'aoi':
        assert counts['point_accounting']['excluded_domain'] > 0, case
    return dict(status='valid', workers=workers, maximum_admitted_workers=max(admitted_workers),
                shared_instance_read_groups=len(groups), comparisons=checks)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    root = args.root
    gdal.UseExceptions()
    gdal.SetConfigOption('GDAL_PAM_ENABLED', 'NO')
    runtime = read_json(root / 'runtime.json')
    with laspy.open(root / 'input/input.copc.laz') as reader:
        count = reader.header.point_count
        extras = list(reader.header.point_format.extra_dimension_names)
    assert extras.index('PredInstance') == 14, extras
    # Fixture tripwire: accidental equal labels would let a shared-reader mixup
    # pass all differential comparisons without exercising independent groups.
    fixture = laspy.read(root / 'input/source.laz')
    label_sets = [set(np.unique(fixture[name])) - {0}
                  for name in ('PredInstance', 'PredInstance_SAT', 'PredInstance_FM')]
    assert all(label_sets), label_sets
    assert all(not left.intersection(right) for i, left in enumerate(label_sets)
               for right in label_sets[i + 1:]), label_sets
    del fixture
    for dimension in ('PredInstance', 'PredInstance_SAT', 'PredInstance_FM'):
        order = read_json(root / f'input/order_{dimension}.json')
        assert order['multiset_equal'] and order['original_point_index_unique'], order
        assert order['original_point_index_restores_expected_order'], order
    resources = {}
    for name in runtime['cases']:
        resource = read_json(root / name / 'resources.json')
        assert resource['exit_code'] == 0, (name, resource)
        assert resource['memory_peak_bytes'] < runtime['memory_bytes'], (name, resource)
        assert resource['memory_limit_bytes'] == runtime['memory_bytes'], (name, resource)
        assert resource['memory_swap_limit_bytes'] == 0, (name, resource)
        assert not resource['memory_events'].get('oom', 0), (name, resource)
        assert not resource['memory_events'].get('oom_kill', 0), (name, resource)
        resource['wall_seconds'] = int((root / name / 'elapsed_seconds.txt').read_text())
        resources[name] = resource
    checks = {f'{scope}-w{workers}': compare_case(root, scope, workers, count)
              for scope in ('full', 'aoi') for workers in (1, 2, 4)}
    report = dict(status='valid', point_count=count, resources=resources, cases=checks,
                  source_sha256=hashlib.sha256((root / 'input/input.copc.laz').read_bytes()).hexdigest(),
                  comparison='Exact CSV bytes, complete GeoJSON and raster cells/metadata; no numeric tolerance',
                  runtime=runtime)
    (root / 'comparison.json').write_text(json.dumps(report, indent=2) + '\n')
    print(f'{len(checks)} parallel production comparisons passed; evidence: {root / "comparison.json"}')


if __name__ == '__main__':
    main()
