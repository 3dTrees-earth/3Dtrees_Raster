"""Compare synthetic production products, allowing only the new output prefix."""
import filecmp
import json
from pathlib import Path
import sys

root = Path(sys.argv[1])
checked = 0
for source in (root / "oracle").iterdir():
    if source.suffix not in (".csv", ".geojson", ".tif") or "performance" in source.name:
        continue
    for case in ("full", "aoi"):
        # Gap filling depends on the footprint; only the AOI case shares the oracle's.
        if case == "full" and source.name.endswith("_chm_filled.tif"):
            continue
        target = root / case / source.name.replace("150_", "foreststructure_", 1)
        if source.suffix == ".geojson":
            expected = json.loads(source.read_text())
            expected["name"] = source.stem.replace("150_", "foreststructure_", 1)
            assert json.loads(target.read_text()) == expected, target
        else:
            assert filecmp.cmp(source, target, shallow=False), target
        checked += 1
print(f"{checked} scientific comparisons passed (only GeoJSON layer prefix changes)")
