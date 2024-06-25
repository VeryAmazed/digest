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
            'build/include',
            'pybind',
        get_pybind_include()
    ],
    library_dirs=['build/lib'],
    define_macros = [("PYBIND", None)],
    extra_compile_args=['-std=c++17', '-fPIC']
)

# class MesonBuildExt(build_ext):
#     def run(self):
#         # Check for Meson and Ninja installation
#         try:
#             subprocess.check_output(['meson', '--version'])
#         except FileNotFoundError:
#             raise RuntimeError("Meson must be installed to build the extensions")
        
#         try:
#             subprocess.check_output(['ninja', '--version'])
#         except FileNotFoundError:
#             raise RuntimeError("Ninja must be installed to build the extensions")

#     def build_extension(self):
#         build_temp = os.path.abspath(self.build_temp)
#         # ext_fullpath = self.get_ext_fullpath(ext.name)
#         # ext_dir = os.path.abspath(os.path.dirname(ext_fullpath))
#         ext_dir = os.path.abspath(os.path.join(ROOT_DIR, 'build'))
#         meson_build_dir = os.path.join(build_temp, 'meson_build')

#         # Create build directory if it doesn't exist
#         if not os.path.exists(meson_build_dir):
#             os.makedirs(meson_build_dir)

#         meson_args = [
#             'meson', 'setup', '--prefix', ext_dir,
#             '--buildtype=release', meson_build_dir
#         ]

#         subprocess.check_call(meson_args)
#         subprocess.check_call(['meson', 'install', '-C', meson_build_dir])

setup(
    name = 'Digest',
    version = '0.2',
    python_requires=">=3.8",
    setup_requires=['setuptools', 'pybind11>=2.6.0', 'meson','ninja'],
    install_requires=['pybind11>=2.6.0'],
    packages=find_packages(),
    ext_modules = [digest],
    include_package_data=True,
    # cmdclass={'build_ext': MesonBuildExt},
)
