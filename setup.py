### Adapted from Uncalled setup.py by Sam Kovaka

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
import os
import subprocess
import sys

try:
    from pybind11.setup_helpers import Pybind11Extension, ParallelCompile, naive_recompile
    ParallelCompile("NPY_NUM_BUILD_JOBS", default=4, needs_recompile=naive_recompile).install()
except:
    from setuptools import Extension as Pybind11Extension

ROOT_DIR = os.getcwd()

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __str__(self):
        import pybind11
        from pybind11.setup_helpers import Pybind11Extension
        return pybind11.get_include()

digest = Pybind11Extension(
    "Digest",
    sources = ['pybind/bindings.cpp'],
    include_dirs = [
            'include',
            'pybind',
        get_pybind_include()
    ],
    define_macros = [("PYBIND", None)]
)

setup(
    name = 'Digest',
    version = '0.2',
    python_requires=">=3.8",
    setup_requires=['setuptools', 'pybind11>=2.6.0'],
    install_requires=['pybind11>=2.6.0'],
    packages=find_packages(),
    ext_modules = [digest],
    include_package_data=True,
)