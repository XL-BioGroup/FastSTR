"""Schedule independent cross-window jobs; restore original output order.

No sequence/window/candidate edits. Long payloads begin first; chunksize=1
prevents unrelated jobs from queuing behind a slow job inside one chunk.
The original payload transfer is retained (fork and spawn both supported).
"""
import multiprocessing
from .make_consensus_TRs import calculate_cross_subread_tr
from .memory_policy import create_pool, release_free_heap


def _calculate(index, args):
    result = calculate_cross_subread_tr(*args)
    release_free_heap(force_check=True)
    return index, result


def calculate_cross_tasks(tasks, workers):
    if not tasks:
        return []
    ordered = sorted(enumerate(tasks), key=lambda item: -len(item[1][1]))
    pool = create_pool(multiprocessing.get_context(), workers)
    try:
        results = pool.starmap(_calculate, ordered, chunksize=1)
        pool.close()
        pool.join()
    except BaseException:
        pool.terminate()
        pool.join()
        raise
    restored = [None] * len(tasks)
    for index, result in results:
        restored[index] = result
    return restored
