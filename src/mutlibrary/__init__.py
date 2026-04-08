"""
mutlibrary – demultiplexing and read processing tools.

This package provides:
- core
- models
- services
- jobs
- api
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("mutlibrary")
except PackageNotFoundError:  # pragma: no cover
    __version__ = "0.1.0"

# from .models.demultiplex import DemultiplexConfig  # noqa: F401
# from .models.read import Read  # noqa: F401
# from .core.demultiplex import (  # noqa: F401
#    build_decision_callbacks,
#    decide_on_barcode,
#    count_start_kmers,
# )
# from .jobs.demultiplex import DemultiplexJob  # noqa: F401

# __all__ = [
#    "__version__",
#    "DemultiplexConfig",
#    "Read",
#    "build_decision_callbacks",
#    "decide_on_barcode",
#    "count_start_kmers",
#    "DemultiplexJob",
# ]
