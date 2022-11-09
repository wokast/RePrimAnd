from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension

ext_modules = [
    Pybind11Extension(
        "pyreprimand",
        ["pyreprimand.cpp"],
        libraries=["RePrimAnd"]
    ),
]

setup(ext_modules=ext_modules)
