The full documentation can be found [here](https://wokast.github.io/RePrimAnd/index.html)

# Installation

The main target and devolpment platform is Linux, although the
library code is not platform-specific and should also work on Macs.
Windows and AIX are not supported.

## Requirements

* A C++11 capable compiler (tested with gcc and clang). 
* Meson build system.
* Boost library.
* GNU Scientific Library version >= 2.0 
* HDF5 library (Only C-bindings required, not the C++ API).
* Doxygen (only for documentation)
* Sphinx with Breathe and bibtex extensions (only for documentation)
* Python matplotlib (only for benchmark plots)
* Python pybind11 package >= 2.6.0 (only for Python bindings)

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

The Python bindings are automatically build together with the library
itself, provided that a Python3 installation is detected which contains
the pybind11 package. The latter can be easily installed via pip or 
conda. Make sure to install pybind11 version >=2.6.0 or the build will fail. 
When using a conda virtual environment, also set the meson `--prefix` option
to the root of the environment, even when you only need the python module.

To disable building the Python bindings and remove the 
corresponding dependencies, use the build option

```bash
meson configure -Dbuild_python_api=false
```

If the Python bindings have been build, they are automatically installed 
together with the library. The Python extension module is called 
`pyreprimand`.


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

Currently there is are no thorns that provide an ET-style interface for using
this library, although this might change in future versions. ET Thorns can of
course use the library like any other C++ library.


## Creating Documentation

To just build the documentation, use the target `documentation`.

```bash
ninja documentation
```

The resulting pages can be found in the build directory under
`docsrc/sphinx/index.html`.

The building of the documentation requires sphinx with the breathe 
and sphinxcontrib-bibtex extensions as well as doxygen.
To disable building documentation and remove 
the corresponding dependencies, use the build option

```bash
meson configure -Dbuild_documentation=false
```

If the documentation is build, it is installed automatically when 
installing the library, by default to 
`<prefix>/usr/local/share/doc/libreprimand/index.html`.


## Running Tests

### Unit tests

Please also take a minute to run the unit tests to ensure 
correct compilation

```bash
ninja test
```

### Benchmarks

The repository contains code to map the efficiency and accuracy of
the primitive recovery, producing the plots shown in the 
article. To recreate the data and plots,


```bash
ninja benchplots
ninja accuracyplots
```

The resulting pdf figures are placed in the build directory under
`tests/benchmarks`.

This requires Python+matplotlib. To disable building benchmarking and 
remove the corresponding dependencies, use the build option

```bash
cd <repository>/mbuild
meson configure -Dbuild_benchmarks=false
```

### Visualizing Master Function

In addition, there is code to sample the primitive recovery master
function (the central ingredient of the scheme) for various cases,
as shown in the paper.

```bash
ninja srootdata
```

The resulting data files are placed in the build directory under 
`tests/sample_root/`.




