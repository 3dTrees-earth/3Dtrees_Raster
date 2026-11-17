"""Run the kNN filter's window parts in worker processes without risking the container.

A part starts when three things hold: a worker is free; the reservations of
the running parts and its own fit the budget; and the container's measured
working set plus those reservations stays below its limit with
MEMORY_FREE_SHARE kept free. A reservation is the part's estimate until parts
of real size (estimated at LEARN_FROM_BYTES or more) have finished; from then
on it is the estimate times their largest measured peak share, times
LEARNED_MARGIN, never above the estimate. So the parts run as close to their
real use as measurements allow, and the live check still holds if an estimate
is wrong. A part that fits nothing beside it runs alone.

If the kernel kills a worker all the same (workers set the highest OOM score
in their initializer), the parts that had not finished run again one at a time
with all threads, then the rest continues in parallel on the plain estimates;
a part killed while running alone stops the run with memory_budget_exceeded.
Finished parts are kept. A part's result must not depend on when or beside
what it runs, so every path gives the same output.

The work function, its worker initializer and the memory probe are passed in,
so tests can run this scheduler with fake parts, a fake container or a worker
that dies.
"""
from __future__ import annotations

from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from concurrent.futures.process import BrokenProcessPool
import multiprocessing

from container_memory import working_set_and_limit

LEARNED_MARGIN = 1.5
LEARN_FROM_BYTES = 1024**3
# Share of the container's memory limit kept free when a part starts.
MEMORY_FREE_SHARE = 0.10
# Seconds between admission checks while parts run: memory frees up meanwhile.
POLL_SECONDS = 5


def run_parts(pending, workers, budget, threads, work, initializer=None, initargs=(),
              memory_probe=working_set_and_limit):
    """Run every (estimate, task) of `pending`; return (results, learned share).

    `work(task + (threads,))` runs in a worker process and returns a tuple whose
    last element is the part's measured peak memory in bytes. Results come in
    completion order. `memory_probe()` returns the container's (working set,
    limit) in bytes, the limit None when there is none.
    """
    pending = list(pending)
    results = []
    learned = 1.0
    retry = []  # parts lost with a killed worker: rerun one at a time

    def reservation(estimate):
        return estimate * learned

    def fits(estimate, running):
        reserved = sum(reservation(e) for e, _ in running.values()) + reservation(estimate)
        if reserved > budget:
            return False
        working, cap = memory_probe()
        return cap is None or working + reserved <= cap * (1 - MEMORY_FREE_SHARE)

    def new_pool(size):
        return ProcessPoolExecutor(size, mp_context=multiprocessing.get_context("spawn"),
                                   initializer=initializer, initargs=initargs,
                                   max_tasks_per_child=1)

    def record(estimate, result):
        nonlocal learned
        peak = result[-1]
        if estimate >= LEARN_FROM_BYTES:
            share = min(1.0, LEARNED_MARGIN * peak / estimate)
            learned = share if learned == 1.0 else max(learned, share)
        results.append(result)

    while pending or retry:
        if retry:
            estimate, task = retry.pop(0)
            single = new_pool(1)
            try:
                result = single.submit(work, task + (max(1, threads),)).result()
            except BrokenProcessPool:
                raise SystemExit("RASTER_FLAG memory_budget_exceeded: a kNN window part was killed "
                                 "for memory while running alone") from None
            except RuntimeError as error:
                raise SystemExit(str(error)) from None
            finally:
                single.shutdown(wait=True, cancel_futures=True)
            record(estimate, result)
            continue
        pool = new_pool(workers)
        running = {}
        try:
            while pending or running:
                while pending and len(running) < workers and (not running or fits(pending[0][0], running)):
                    estimate, task = pending.pop(0)
                    # A part that leaves room for few others queries with the idle threads.
                    fit = max(1, min(workers, int(budget // max(reservation(estimate), 1))))
                    future = pool.submit(work, task + (max(1, threads // fit),))
                    running[future] = (estimate, task)
                done, _ = wait(running, timeout=POLL_SECONDS, return_when=FIRST_COMPLETED)
                for future in done:
                    estimate, task = running[future]
                    try:
                        result = future.result()
                    except BrokenProcessPool:
                        # A killed worker (BrokenProcessPool is a RuntimeError): retried below.
                        raise
                    except RuntimeError as error:
                        raise SystemExit(str(error)) from None
                    # Removed only once its result is in: a killed worker's part is retried.
                    del running[future]
                    record(estimate, result)
        except BrokenProcessPool:
            retry = list(running.values())
            print(f"KNN_RETRY parts={len(retry)} reason=worker_killed mode=one_at_a_time", flush=True)
            # The rest goes on in parallel, on the plain estimates until parts finish again.
            learned = 1.0
        finally:
            pool.shutdown(wait=True, cancel_futures=True)
    return results, learned
