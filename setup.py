# Welcome to the OpenBNSL setup.py.
# Environment variables you are probably interested in:
# TODO
#   -Default: TBD
#   -Description: TBD
# 

import os
import subprocess
import sys
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext

class CMakeBuildExtension(build_ext):
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            "-DCMAKE_BUILD_TYPE=Release",
        ]

        build_temp = self.build_temp
        os.makedirs(build_temp, exist_ok=True)

        try:
            subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=build_temp)
            subprocess.check_call(["cmake", "--build", "."], cwd=build_temp)
        except subprocess.CalledProcessError as e:
            print(f"Error during CMake build: {e}")
            raise

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

description = "A Python package for the OpenBNSL core library."
try:
    long_description = open("README.md").read()
except FileNotFoundError:
    long_description = description

setup(
    name="openbnsllib",
    version="0.1.0",
    description=description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    ext_modules=[CMakeExtension("openbnsllib")],
    cmdclass={"build_ext": CMakeBuildExtension},
    packages=find_packages(include=["helpers", "helpers.*"]),
    install_requires=[], # TODO: Add dependencies
    extras_require={},
    url="https://github.com/hal-lab-u-tokyo/openbnsl",
    author="Ryota Miyagi",
    author_email="rmiyagi@hal.ipc.i.u-tokyo.ac.jp",
    python_requires=">=3.7",
    license = "MIT",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "Topic :: Software Development",
        "Topic :: Software Development :: Libraries",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Programming Language :: C++",
        "Programming Language :: Python :: 3",
    ],
    keywords="Bayesian Network, Structure Learning, Machine Learning",

)
