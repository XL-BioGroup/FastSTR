"""Build optional native accelerators; Python implementations remain available."""
import os

from Cython.Build import cythonize
from setuptools import Extension, setup


mode = os.environ.get("FASTSTR_BUILD_EXTENSIONS", "auto")
if mode not in ("auto", "0", "1"):
    raise ValueError("FASTSTR_BUILD_EXTENSIONS must be auto, 0, or 1")

extensions = []
if mode != "0":
    extensions = cythonize(
        [Extension("faststr." + name, ["faststr/" + name + ".pyx"],
                   optional=(mode == "auto"))
         for name in ("_anchors", "_clusters", "_motif_stats", "_rpta", "_dbscan")],
        compiler_directives={"language_level": 3},
        build_dir="build/cython",
    )
    # Cython recreates Extension objects and can drop setuptools' optional flag.
    for extension in extensions:
        extension.optional = mode == "auto"

setup(ext_modules=extensions)
