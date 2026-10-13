"""Memory housekeeping only: preserve all detector parameters and results.

Freeze pre-fork Python objects in workers to avoid cyclic-GC COW writes.
The parent immediately restores its original freeze state after pool creation.
Occasionally return already-free glibc heap pages; never discard live objects.
"""
import ctypes
import gc
import os

_ENABLED = os.environ.get('FASTSTR_MEMORY_POLICY', 'on') != 'off'
_calls = 0
_libc = None
_trim = None
_info = None


class _Mallinfo(ctypes.Structure):
    _fields_ = [(name, ctypes.c_size_t) for name in
                ('arena', 'ordblks', 'smblks', 'hblks', 'hblkhd', 'usmblks',
                 'fsmblks', 'uordblks', 'fordblks', 'keepcost')]


if _ENABLED:
    try:
        _libc = ctypes.CDLL(None)
        _trim = _libc.malloc_trim
        _trim.argtypes = [ctypes.c_size_t]
        _trim.restype = ctypes.c_int
        _info = _libc.mallinfo2
        _info.argtypes = []
        _info.restype = _Mallinfo
    except (AttributeError, OSError):
        _trim = _info = None


def create_pool(context, workers):
    """Do not alter user-owned frozen objects or automatic-GC enabled state."""
    owned_freeze = (_ENABLED and context.get_start_method() == 'fork'
                    and hasattr(gc, 'freeze') and gc.get_freeze_count() == 0)
    if owned_freeze:
        gc.freeze()
    try:
        return context.Pool(workers)
    finally:
        if owned_freeze:
            gc.unfreeze()


def release_free_heap(force_check=False):
    """Housekeeping cadence/bytes threshold, NOT STR detection thresholds."""
    global _calls
    if not _ENABLED or _trim is None:
        return
    _calls += 1
    if (force_check or _calls % 128 == 0) and _info().fordblks >= 32 * 1024 * 1024:
        _trim(0)
