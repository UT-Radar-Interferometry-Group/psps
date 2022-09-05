"""
Based on https://github.com/pybind/python_example
"""
from setuptools import setup
import sys
import sysconfig
import os
from glob import glob

curpath = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(curpath, "extern", "pybind11"))
from pybind11.setup_helpers import Pybind11Extension, build_ext


__version__ = "0.0.1"

# The main interface is through Pybind11Extension.
# * You can add cxx_std=11/14/17, and then build_ext can be removed.
# * You can set include_pybind11=false to add the include directory yourself,
#   say from a submodule.
#
# Note:
#   Sort input source files if you glob sources to ensure bit-for-bit
#   reproducible builds (https://github.com/pybind/python_example/pull/53)


ext = Pybind11Extension(
    "psps.cppsimilarity",
    glob("src/cpp/*.cpp"),
    # passing in the version to the compiled code
    define_macros=[("VERSION_INFO", __version__)],
)

# start with python level flags
extra_compile_args = sysconfig.get_config_var("CXXFLAGS")
extra_compile_args = extra_compile_args.split() if extra_compile_args else []
# Note: passing openmp as "compile_arg" for the separate pieces, and a 'link_arg'
# for the final linking that setuptools does to make the python .so object
# setuptools also seems to pick which std c++ to use
extra_compile_args += ["-Wall", "-Wextra", "-fopenmp"]  # -std=c++11
extra_link_args = ["-fopenmp"]
# Mac doesnt seem to allow -fopenmp for me
if sys.platform == "darwin":
    extra_link_args.pop(-1)
    extra_compile_args.pop(-1)
    # IDK how to make this work on Mac yet... does not work:
    # https://stackoverflow.com/questions/39979836/using-openmp-with-c11-on-mac-os

ext.extra_compile_args += extra_compile_args
ext.extra_link_args += extra_link_args

setup(
    name="psps",
    version=__version__,
    author="Ke Wang",
    author_email="kewang@utexas.edu",
    url="https://github.com/aik9508/PSInSAR",
    description="Persistent Scatterer selection based on phase similarity",
    long_description="",
    ext_modules=[ext],
    cmdclass={"build_ext": build_ext},
    install_requires=[
        "matplotlib",
        "numpy",
        "scipy",
        "notebook"
    ],
    zip_safe=False,
)
