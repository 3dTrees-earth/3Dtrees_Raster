"""Initial Galaxy allocation from the merged COPC's recorded size in bytes.

This dispatch policy does not open point clouds or estimate local point density.
The analysis runtime independently admits spatial tasks within the allocation.
"""

import argparse
from dataclasses import asdict, dataclass
import json


POLICY_VERSION = "merged-copc-size-2026-10-13"
MIB = 1024 ** 2


@dataclass(frozen=True)
class ResourceRequest:
    policy: str
    input_bytes: int
    tier: str
    cpus: int
    memory_mib: int

    @property
    def memory_gib(self):
        """TPV uses GiB; preserve the integer MiB request exactly."""
        return self.memory_mib / 1024


def request_resources(input_bytes):
    """Return a monotone request capped at 10 CPUs and 100 decimal GB.

    input_bytes must describe the single merged, labelled COPC supplied to this
    job, excluding AOIs, original uploads, and decompression multipliers. Missing
    or invalid metadata fails explicitly rather than requesting a small default.
    File size is an initial allocation heuristic, not a peak-memory guarantee.
    """
    if isinstance(input_bytes, bool) or not isinstance(input_bytes, int) or input_bytes <= 0:
        raise ValueError("merged COPC input_bytes must be a positive integer")
    if input_bytes <= 1_000_000_000:
        tier, cpus, memory_gb = "up_to_1gb", 8, 30
    elif input_bytes <= 5_000_000_000:
        tier, cpus, memory_gb = "up_to_5gb", 10, 30
    elif input_bytes <= 10_000_000_000:
        tier, cpus, memory_gb = "up_to_10gb", 10, 60
    else:
        tier, cpus, memory_gb = "over_10gb", 10, 100
    # Round down so a scheduler's integer MiB request never exceeds the GB cap.
    return ResourceRequest(POLICY_VERSION, input_bytes, tier, cpus,
                           memory_gb * 1_000_000_000 // MIB)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-bytes", required=True, type=int,
                        help="Actual merged COPC size, from Galaxy dataset metadata.")
    args = parser.parse_args()
    try:
        request = request_resources(args.input_bytes)
    except ValueError as error:
        parser.error(str(error))
    print(json.dumps({**asdict(request), "memory_gib": request.memory_gib}, sort_keys=True))


if __name__ == "__main__":
    main()
