#!/usr/bin/env python3
"""Conservative spatial point counts without reading compressed point records.

COPC 1.0 hierarchy entries partition the full source points across *all* depths,
not just leaves. Sum every intersecting data node, including ancestors. Closed
XY intersections deliberately overestimate reads at edges and include all Z.
The caller supplies the actual buffered read windows, in source coordinates.

Only fixed-size LAS/COPC headers and hierarchy entries are read. The fixed
metadata/window limits bound work and auxiliary memory independently of source
point count. Invalid COPC metadata fails closed, never falling back silently.
Format reference: https://copc.io/ (COPC 1.0).
"""

import argparse
import json
import math
import os
from pathlib import Path
import struct
import tempfile


ENTRY = struct.Struct("<iiiiQii")
MAX_METADATA_BYTES = 64 * 1024 * 1024
MAX_WINDOWS = 100_000
MAX_WINDOWS_BYTES = 32 * 1024 * 1024


class MetadataError(ValueError):
    """Input cannot provide a trustworthy conservative count."""


def _intersects(a, b):
    return a[0] <= b[2] and b[0] <= a[2] and a[1] <= b[3] and b[1] <= a[3]


def _windows(value):
    if not isinstance(value, list) or len(value) > MAX_WINDOWS:
        raise MetadataError(f"windows must be a list of at most {MAX_WINDOWS} bounds")
    result = []
    for index, item in enumerate(value):
        if not isinstance(item, (list, tuple)) or len(item) != 4:
            raise MetadataError(f"window {index} must be [xmin, ymin, xmax, ymax]")
        if any(isinstance(v, bool) or not isinstance(v, (int, float)) for v in item):
            raise MetadataError(f"window {index} contains a nonnumeric bound")
        try:
            item = tuple(float(v) for v in item)
        except OverflowError as exc:
            raise MetadataError(f"window {index} exceeds supported coordinate range") from exc
        if not all(math.isfinite(v) for v in item) or item[0] > item[2] or item[1] > item[3]:
            raise MetadataError(f"window {index} has nonfinite or reversed bounds")
        result.append(item)
    return result


class _WindowIndex:
    """Small static bounding-volume tree; no nodes-by-windows dense array."""

    def __init__(self, windows, query_bounds):
        self.windows = windows
        self.counts = [0] * len(windows)
        indices = [i for i, window in enumerate(windows) if _intersects(window, query_bounds)]
        self.root = self._build(indices) if indices else None

    def _build(self, indices):
        bounds = (min(self.windows[i][0] for i in indices),
                  min(self.windows[i][1] for i in indices),
                  max(self.windows[i][2] for i in indices),
                  max(self.windows[i][3] for i in indices))
        if len(indices) <= 16:
            return bounds, indices, None, None
        axis = 0 if bounds[2] - bounds[0] >= bounds[3] - bounds[1] else 1
        indices.sort(key=lambda i: self.windows[i][axis] / 2 + self.windows[i][axis + 2] / 2)
        middle = len(indices) // 2
        return bounds, None, self._build(indices[:middle]), self._build(indices[middle:])

    def add(self, bounds, count):
        stack = [self.root] if self.root else []
        while stack:
            node = stack.pop()
            if not _intersects(bounds, node[0]):
                continue
            if node[1] is None:
                stack.extend((node[2], node[3]))
            else:
                for index in node[1]:
                    if _intersects(bounds, self.windows[index]):
                        self.counts[index] += count


class _MetadataReader:
    def __init__(self, stream, limit):
        self.stream = stream
        self.size = os.fstat(stream.fileno()).st_size
        self.limit = limit
        self.bytes_read = 0

    def range(self, offset, size, name):
        if offset < 0 or size <= 0 or offset + size > self.size:
            raise MetadataError(f"{name} range is empty or outside file bounds")

    def read(self, offset, size):
        self.range(offset, size, "metadata")
        if self.bytes_read + size > self.limit:
            raise MetadataError(f"hierarchy exceeds {self.limit} byte metadata limit")
        self.stream.seek(offset)
        data = self.stream.read(size)
        self.bytes_read += len(data)
        if len(data) != size:
            raise MetadataError("truncated metadata read")
        return data


def _header(reader):
    raw = reader.read(0, 227)
    if raw[:4] != b"LASF" or raw[24] != 1 or raw[25] not in (2, 3, 4):
        raise MetadataError("expected a LAS 1.2, 1.3 or 1.4 header")
    minor = raw[25]
    required_size = {2: 227, 3: 235, 4: 375}[minor]
    if required_size > len(raw):
        raw += reader.read(len(raw), required_size - len(raw))
    header_size, points_offset, vlr_count = struct.unpack_from("<HII", raw, 94)
    if header_size < required_size or not header_size <= points_offset <= reader.size:
        raise MetadataError("invalid LAS header size or point data offset")
    legacy_count = struct.unpack_from("<I", raw, 107)[0]
    extended_count = struct.unpack_from("<Q", raw, 247)[0] if minor == 4 else 0
    # Plain LAS 1.4 readers can use the legacy count when the extended field is
    # zero. Conflicting nonzero fields must never lower a memory reservation.
    point_count = max(legacy_count, extended_count)
    max_x, min_x, max_y, min_y, max_z, min_z = struct.unpack_from("<6d", raw, 179)
    mins, maxs = (min_x, min_y, min_z), (max_x, max_y, max_z)
    if not all(math.isfinite(v) for v in mins + maxs) or any(a > b for a, b in zip(mins, maxs)):
        raise MetadataError("invalid LAS coordinate bounds")
    scales = struct.unpack_from("<3d", raw, 131)
    if not all(math.isfinite(v) and v > 0 for v in scales):
        raise MetadataError("invalid LAS coordinate scales")
    # Writers may retain pre-quantization extrema in the public header. A
    # stored point can then lie just beyond those bounds, including exactly on
    # a catalog boundary. Use one coordinate quantum and outward rounding for
    # conservative pruning; keep the raw header bounds as provenance only.
    query_bounds = (
        math.nextafter(min_x - scales[0], -math.inf),
        math.nextafter(min_y - scales[1], -math.inf),
        math.nextafter(max_x + scales[0], math.inf),
        math.nextafter(max_y + scales[1], math.inf),
    )
    if not all(math.isfinite(v) for v in query_bounds):
        raise MetadataError("LAS query bounds exceed supported coordinate range")
    info = None
    if vlr_count:
        if header_size + 54 > points_offset:
            raise MetadataError("first VLR overlaps point records")
        first = reader.read(header_size, 54)
        user_id = first[2:18].split(b"\0", 1)[0]
        record_id, length = struct.unpack_from("<HH", first, 18)
        if header_size + 54 + length > points_offset:
            raise MetadataError("first VLR extends into point records")
        if user_id == b"copc" and record_id == 1:
            if minor != 4 or header_size != 375 or length != 160 or (raw[104] & 0x3f) not in (6, 7, 8) or not raw[104] & 0x80:
                raise MetadataError("invalid COPC 1.0 header, info VLR or point format")
            # COPC always uses the extended count; the hierarchy must reconcile
            # against that field, including when it is zero. Existing canonical
            # writers can retain a nonzero deprecated legacy counter; ignore it.
            point_count = extended_count
            info = reader.read(header_size + 54, 160)
    return point_count, mins, maxs, points_offset, info, query_bounds


def _node_bounds(key, center, halfsize):
    level, x, y, _ = key
    width = math.ldexp(halfsize, 1 - level)
    if width == 0:
        raise MetadataError("hierarchy depth underflows spatial bounds")
    # Round outwards to avoid losing nodes at representational boundaries.
    result = []
    for axis, coordinate in ((0, x), (1, y)):
        low = center[axis] - halfsize + width * coordinate
        high = center[axis] - halfsize + width * (coordinate + 1)
        tolerance = 4 * max(math.ulp(center[axis]), math.ulp(halfsize), math.ulp(low), math.ulp(high))
        result.append((low - tolerance, high + tolerance))
    return result[0][0], result[1][0], result[0][1], result[1][1]


def estimate(input_path, windows, *, allow_header_fallback=False, max_metadata_bytes=MAX_METADATA_BYTES):
    """Return upper bounds in input-window order, or raise MetadataError.

    The optional non-COPC fallback returns the complete header count for every
    intersecting window; it never estimates point density from compressed size.
    The source must remain immutable for the duration of this metadata scan.
    """
    windows = _windows(windows)
    with open(input_path, "rb") as stream:
        reader = _MetadataReader(stream, max_metadata_bytes)
        point_count, mins, maxs, points_offset, info, query_bounds = _header(reader)
        header_bounds = (mins[0], mins[1], maxs[0], maxs[1])
        result = dict(schema_version=1, point_count=point_count, windows=windows,
                      header_bounds=list(header_bounds), query_bounds=list(query_bounds),
                      metadata_bytes_read=reader.bytes_read)
        if info is None:
            if not allow_header_fallback:
                raise MetadataError("COPC info VLR not found; header fallback requires explicit opt-in")
            result.update(mode="header_fallback", point_count_upper_bounds=[
                point_count if _intersects(w, query_bounds) else 0 for w in windows])
            return result
        center = struct.unpack_from("<3d", info)
        halfsize, spacing = struct.unpack_from("<dd", info, 24)
        if not all(math.isfinite(v) for v in center + (halfsize, spacing)) or halfsize <= 0 or spacing <= 0 or any(info[72:]):
            raise MetadataError("invalid COPC center, size, spacing or reserved fields")
        if not math.isfinite(halfsize * 2):
            raise MetadataError("COPC root width exceeds supported coordinate range")
        for axis in range(3):
            low, high = center[axis] - halfsize, center[axis] + halfsize
            if not math.isfinite(low) or not math.isfinite(high):
                raise MetadataError("COPC root bounds overflow")
            tolerance = 4 * max(math.ulp(low), math.ulp(high), math.ulp(halfsize))
            if mins[axis] < low - tolerance or maxs[axis] > high + tolerance:
                raise MetadataError("COPC root cube does not enclose LAS header bounds")
        root_offset, root_size = struct.unpack_from("<QQ", info, 40)
        pending = [(root_offset, root_size)]
        pages, keys, ranges = set(), set(), []
        total, entries = 0, 0
        index = _WindowIndex(windows, query_bounds)
        while pending:
            offset, size = pending.pop()
            if offset in pages:
                raise MetadataError("duplicate hierarchy page reference or page cycle")
            reader.range(offset, size, "hierarchy page")
            if offset < points_offset or size % ENTRY.size:
                raise MetadataError("invalid hierarchy page offset or entry alignment")
            if reader.bytes_read + size > reader.limit:
                raise MetadataError(f"hierarchy exceeds {reader.limit} byte metadata limit")
            pages.add(offset)
            ranges.append((offset, offset + size))
            # A page's declared size never becomes a single potentially huge allocation.
            for start in range(offset, offset + size, 64 * 1024):
                data = reader.read(start, min(64 * 1024, offset + size - start))
                for level, x, y, z, child_offset, byte_size, count in ENTRY.iter_unpack(data):
                    entries += 1
                    # int32 coordinates limit useful keys to 31 bits. Deeper
                    # keys remain legal but cap depth before exponentiation.
                    if not 0 <= level <= 1023 or min(x, y, z) < 0 or max(x, y, z) >= 2 ** level:
                        raise MetadataError("invalid hierarchy voxel key")
                    if count < -1:
                        raise MetadataError("invalid negative hierarchy point count")
                    key = (level, x, y, z)
                    if count == -1:
                        reader.range(child_offset, byte_size, "child hierarchy page")
                        pending.append((child_offset, byte_size))
                        continue
                    packed_key = struct.pack("<iiii", *key)
                    if packed_key in keys:
                        raise MetadataError("duplicate hierarchy data key")
                    keys.add(packed_key)
                    if count == 0:
                        if child_offset != 0 or byte_size != 0:
                            raise MetadataError("empty hierarchy node has a nonempty data range")
                        continue
                    reader.range(child_offset, byte_size, "point chunk")
                    if child_offset < points_offset:
                        raise MetadataError("point chunk overlaps LAS headers")
                    ranges.append((child_offset, child_offset + byte_size))
                    total += count
                    if total > point_count:
                        raise MetadataError("hierarchy point count exceeds LAS header count")
                    index.add(_node_bounds(key, center, halfsize), count)
        if total != point_count:
            raise MetadataError(f"hierarchy point count {total} differs from LAS header {point_count}")
        ranges.sort()
        if any(a[1] > b[0] for a, b in zip(ranges, ranges[1:])):
            raise MetadataError("overlapping hierarchy pages or point chunks")
        result.update(mode="copc_hierarchy", hierarchy_point_count=total,
                      hierarchy_node_count=len(keys), hierarchy_entry_count=entries,
                      hierarchy_page_count=len(pages), point_count_upper_bounds=index.counts,
                      metadata_bytes_read=reader.bytes_read)
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--windows", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--allow-header-fallback", action="store_true")
    args = parser.parse_args()
    if args.output.resolve() in (args.input.resolve(), args.windows.resolve()):
        parser.error("output must differ from input and windows paths")
    try:
        with args.windows.open("rb") as source:
            raw = source.read(MAX_WINDOWS_BYTES + 1)
        if len(raw) > MAX_WINDOWS_BYTES:
            raise MetadataError("windows JSON exceeds size limit")
        result = estimate(args.input, json.loads(raw), allow_header_fallback=args.allow_header_fallback)
        temporary = None
        try:
            with tempfile.NamedTemporaryFile("w", dir=args.output.parent, delete=False) as target:
                temporary = target.name
                json.dump(result, target, allow_nan=False)
                target.write("\n")
            os.replace(temporary, args.output)
        finally:
            if temporary and os.path.exists(temporary):
                os.unlink(temporary)
    except (OSError, ValueError, struct.error) as exc:
        parser.exit(1, f"COPC task count failed: {exc}\n")


if __name__ == "__main__":
    main()
