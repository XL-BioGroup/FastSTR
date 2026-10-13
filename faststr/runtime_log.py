"""Small per-run summaries; no sequence scans or background resource sampling."""
from contextlib import contextmanager
import logging
import os
import time


def _optional_probe(probe):
    try:
        return probe()
    except (AttributeError, OSError, NotImplementedError):
        return "unavailable"


def log_input_file(path):
    logging.info("Input file: path=%s size_bytes=%s", os.path.abspath(path),
                 _optional_probe(lambda: os.path.getsize(path)))


def _log_load(label):
    # OS load averages, not FastSTR CPU utilization. Unsupported on Windows.
    logging.info("System load average %s (1/5/15 min): %s", label,
                 _optional_probe(lambda: os.getloadavg()))


class RunStatus:
    def __init__(self):
        self.stage_name = "initialization"

    def stage(self, name):
        self.stage_name = name
        logging.info("Stage started: %s", name)


@contextmanager
def run_logging(path, workers):
    # Use a dedicated handler even when an API caller configured root logging.
    # Restore the caller's configuration, including on failure or repeat calls.
    handler = logging.FileHandler(path, mode="w", encoding="utf-8")
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
    logger = logging.getLogger()
    previous_level = logger.level
    logger.addHandler(handler)
    logger.setLevel(min(previous_level, logging.INFO))
    started = time.perf_counter()
    status = RunStatus()
    try:
        logging.info("Program started")
        logging.info(
            "Run resources: logical_cpus=%s cpu_affinity_count=%s requested_workers=%s",
            _optional_probe(os.cpu_count),
            _optional_probe(lambda: len(os.sched_getaffinity(0))), workers)
        _log_load("at start")
        yield status
    except (Exception, KeyboardInterrupt):
        logging.exception("Program failed during stage: %s", status.stage_name)
        raise
    else:
        logging.info("Program completed")
    finally:
        try:
            logging.info("Run elapsed wall time: %.6f seconds (log setup to completion/failure; excludes Python imports)",
                         time.perf_counter() - started)
            _log_load("at end")
        finally:
            logger.removeHandler(handler)
            handler.close()
            logger.setLevel(previous_level)
