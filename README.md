# Overview

RePrimAnd is a support library for numerical simulations of general 
relativistic magnetohydrodynamics and other neutron star related
problems, written in C++ (a Python interface is also available)

The full documentation can be found [here](https://wokast.github.io/RePrimAnd/index.html)

## Installation

The main target and devolpment platform is Linux, although the
library code is not platform-specific and should also work on Macs.
Windows and AIX are not supported. For use in HPC, users need to build 
the library from source. For use in postprocessing, the Python intferface
can also be installed via pip (binary wheel of the library available for 
Linux). Conda packages for the C++ library and Python intface are planned 
but not available yet.

## Requirements

* A C++11 capable compiler (tested with gcc and clang). 
* Meson build system.
* Boost library.
* GNU Scientific Library version >= 2.0 
* HDF5 library (Only C-bindings required, not the C++ API).
* Doxygen (only for documentation)
* Sphinx with Breathe and bibtex extensions (only for documentation)
* Python matplotlib (only for benchmark plots)

## Building from Source

The library is build using the modern
[Meson](https://mesonbuild.com>)
build system, which is 
[available](https://mesonbuild.com/Getting-meson.html)
as standard package in most distributions and can also be installed 
via Python's pip.

To build the library, 

```bash
cd <repository>
meson mbuild --buildtype=release --prefix=<custom install location>
cd mbuild
ninja
```

This will compile with optimization and without debug symbols. Other
possibilities are `--buildtype=debug` and `--buildtype=debugoptimized`
To use a different compiler, e.g. clang, prefix the meson command
with `CC=clang CXX=clang++`.
See [here](https://mesonbuild.com/Running-Meson.html) for general 
Meson usage.


## Installing

To install the library, use

```bash
cd mbuild
ninja install
```

This will install in a user-defined location if the `--prefix` option
was given during the build setup, otherwise systemwide. 


## Using the Library

The various header files needed for using the library are installed 
in a subdirectory `reprimand`. The executables have to be linked
with library `RePrimAnd`.

A minimal example can be found in `example/minimal.cc`. Assuming
the library is installed where the compiler (below we use gcc as 
example) can find it, compilation should be straightforward:

```bash
g++ -lRePrimAnd --std=c++11 minimal.cc
```

## Python Bindings

If only the Python interface is required on a Linux platform, 
it is easiest to install from pypi

```bash
pip install reprimand
```

Otherwise, one first has to buid and install the C++ library as shown above.

To build and install the Python interface, do

```bash
cd bindings/python
pip install .
```

This will also pip install packages numpy and pybind11>=2.6.0. 
When using conda environment, it may be better to install those 
first using conda.


The Python extension module is called  `pyreprimand`.


## Einstein Toolkit Support

RePrimAnd does provide a thorn that builds the library within
an EinsteinToolkit (ET) environment, using the ExternalLibraries mechanism. The
thorn can be found in the folder `ET_interface/thorns/RePrimAnd/`. The thorn
depends on the HDF5, GSL, and BOOST ExternalLibraries thorns. Building it via
the ET build system does not require meseon. Note this only builds the library,
but not the tests and Python bindings. 

The thorn is also part of the official ET framework. The version in the master 
brach will typically by ahead of the ET version, but is not guaranteed to be 
stable or compatible.

There are two experimental (and largely undocumented for now) thorns 
that aim to simplifying the usage.
`RePrimAnd_Global_EOS` provides a centralized selection of a global thermal 
EOS (e.g. for evolution and analysis) and a global barotropic EOS (e.g. for 
initial data).
`RePrimAnd_EOS_Omni_API` provides the most important subset of the `EOS_Omni` thorn 
interface, forwarding EOS calls to the reprimand EOS set by `RePrimAnd_Global_EOS`.
It is intended for transitioning existing code to the new interface.


## Creating Documentation

Building the documentation is deactivated by default. 
To build it, specify

```bash
meson configure -Dbuild_documentation=true
```

and then build the target `documentation`.

```bash
ninja documentation
```

The resulting pages can be found in the build directory under
`docsrc/sphinx/index.html`.

The building of the documentation requires sphinx with the breathe 
and sphinxcontrib-bibtex extensions as well as doxygen.

If the documentation is build, it is installed automatically when 
installing the library, by default to 
`<prefix>/usr/local/share/doc/libreprimand/index.html`.


## Running Tests

### Unit tests


Building the unit tests is deactivated by default. 
To build them, specify

```bash
meson configure -Dbuild_tests=true
```
and then

```bash
ninja test
```

Please report errors on the issue tracker.

### Benchmarks

The repository contains code to map the efficiency and accuracy of
the primitive recovery, producing the plots shown in the 
article. To recreate the data and plots,


```bash
meson configure -Dbuild_benchmarks=true
ninja benchplots
ninja accuracyplots
```

The resulting pdf figures are placed in the build directory under
`tests/benchmarks`.

This requires Python+matplotlib. 

### Visualizing Master Function

In addition, there is code to sample the primitive recovery master
function (the central ingredient of the con2prim scheme) for various cases,
as shown in the paper.

```bash
meson configure -Dbuild_benchmarks=true
ninja srootdata
```

The resulting data files are placed in the build directory under 
`tests/sample_root/`.




