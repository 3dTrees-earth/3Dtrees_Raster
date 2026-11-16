#!/usr/bin/env python3
"""knn_parts.run_parts with fake parts: admission by the memory probe, a killed
worker's retry, a part killed while alone, a failing part, and learning."""
import os
import signal
import sys
import tempfile
import time

sys.path.insert(0, os.environ.get("RASTER_SOURCE", "/opt/raster"))
import knn_parts  # noqa: E402

GIB = 1024**3


def part(task):
    """Fake part: (name, start, end, peak); sleeps so overlaps can be seen."""
    name, seconds, peak, mode, marker, threads = task
    if mode == "kill_once":
        try:
            os.close(os.open(marker, os.O_CREAT | os.O_EXCL | os.O_WRONLY))
            os.kill(os.getpid(), signal.SIGKILL)
        except FileExistsError:
            pass
    elif mode == "kill_always":
        os.kill(os.getpid(), signal.SIGKILL)
    elif mode == "fail":
        raise RuntimeError("RASTER_FLAG memory_budget_exceeded: fake part too dense")
    start = time.time()
    time.sleep(seconds)
    return name, start, time.time(), peak


def overlap(results):
    events = sorted([(r[1], 1) for r in results] + [(r[2], -1) for r in results])
    level = most = 0
    for _, step in events:
        level += step
        most = max(most, level)
    return most


def parts(count, estimate=GIB, peak=GIB // 10, mode="", marker=""):
    return [(estimate, (f"p{i}", 0.4, peak, mode if i == 0 else "", marker)) for i in range(count)]


def main():
    free = lambda: (0, None)  # noqa: E731
    results, _ = knn_parts.run_parts(parts(6), 3, 10 * GIB, 3, part, memory_probe=free)
    assert sorted(r[0] for r in results) == [f"p{i}" for i in range(6)], results
    assert overlap(results) == 3, overlap(results)

    # A container near its limit admits one part at a time, never none.
    full = lambda: (9 * GIB, 10 * GIB)  # noqa: E731
    results, _ = knn_parts.run_parts(parts(4), 4, 10 * GIB, 4, part, memory_probe=full)
    assert len(results) == 4 and overlap(results) == 1, overlap(results)

    # Reservations that exceed the budget together run one by one as well: parts that really peak at their
    # estimate keep their reservation.
    results, _ = knn_parts.run_parts(parts(3, estimate=6 * GIB, peak=6 * GIB), 3, 10 * GIB, 3, part,
                                     memory_probe=free)
    assert len(results) == 3 and overlap(results) == 1, overlap(results)

    # Parts that peak far below their estimate teach a small reservation: after the first one the other two
    # run together.
    results, learned = knn_parts.run_parts(parts(3, estimate=6 * GIB), 3, 10 * GIB, 3, part, memory_probe=free)
    assert len(results) == 3 and overlap(results) == 2 and learned < 0.1, (overlap(results), learned)

    # Finished parts of real size teach the reservation share: peak 0.1 of the estimate * 1.5.
    _, learned = knn_parts.run_parts(parts(2), 1, 10 * GIB, 1, part, memory_probe=free)
    assert abs(learned - 0.15) < 1e-9, learned

    with tempfile.TemporaryDirectory() as scratch:
        marker = os.path.join(scratch, "killed")
        results, _ = knn_parts.run_parts(parts(5, mode="kill_once", marker=marker), 3, 10 * GIB, 3, part,
                                         memory_probe=free)
        assert os.path.exists(marker)
        assert sorted(r[0] for r in results) == [f"p{i}" for i in range(5)], results

    for mode, flag in (("kill_always", "killed for memory while running alone"), ("fail", "fake part too dense")):
        try:
            knn_parts.run_parts(parts(3, mode=mode), 2, 10 * GIB, 2, part, memory_probe=free)
        except SystemExit as error:
            assert "RASTER_FLAG memory_budget_exceeded" in str(error) and flag in str(error), error
        else:
            raise AssertionError(f"{mode} did not stop the run")
    print("knn part scheduling checks passed")


if __name__ == "__main__":
    main()
