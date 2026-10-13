"""Same ordered window tasks; fork workers inherit immutable Seq inputs.

Only integer indexes cross the task queue under fork. Non-fork contexts retain
the old sequence-payload path. Pool.map returns original window order; finer
fork-index scheduling changes task granularity, never the detection windows.
"""
import logging
import multiprocessing
import os
import threading
from .get_subread_trs import get_subread_trs
from .memory_policy import create_pool, release_free_heap

_READS = None
_PARAMS = None
_LOCK = threading.Lock()


def _detect_index(index):
    result = get_subread_trs(_READS[index], *_PARAMS)
    release_free_heap()
    return result


def detect_windows(sub_reads, parameters, workers, transport=None):
    global _READS, _PARAMS
    requested = transport or os.environ.get('FASTSTR_TRANSPORT', 'fork_index')
    if requested not in ('fork_index', 'legacy'):
        raise ValueError('FASTSTR_TRANSPORT must be fork_index or legacy')
    context = multiprocessing.get_context()
    use_index = requested == 'fork_index' and context.get_start_method() == 'fork'
    groups = int(os.environ.get('FASTSTR_CHUNKS_PER_WORKER', '16' if use_index else '4'))
    if groups < 1:
        raise ValueError('FASTSTR_CHUNKS_PER_WORKER must be positive')
    chunksize = max(1, (len(sub_reads) + workers * groups - 1) // (workers * groups))
    logging.info('Window transport: %s; context: %s',
                 'fork_index' if use_index else 'legacy', context.get_start_method())
    logging.info('Window task chunksize: %s; groups per worker: %s', chunksize, groups)
    # Lock only protects concurrent parent calls; workers never acquire it.
    with _LOCK:
        if use_index:
            _READS, _PARAMS = sub_reads, tuple(parameters)
        pool = None
        try:
            pool = create_pool(context, workers)
            if use_index:
                result = pool.map(_detect_index, range(len(sub_reads)), chunksize=chunksize)
            else:
                tasks = [(sequence, *parameters) for sequence in sub_reads]
                result = pool.starmap(get_subread_trs, tasks, chunksize=chunksize)
            pool.close()
            pool.join()
            return result
        except BaseException:
            if pool is not None:
                pool.terminate()
                pool.join()
            raise
        finally:
            _READS = _PARAMS = None
