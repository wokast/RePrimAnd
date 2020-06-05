Installation
============

The main target and devolpment platform is Linux, although the 
library code is not platform-specific and should also work on Macs.
Windows and AIX are not supported.

Requirements
------------

* A C++11 capable compiler (tested with gcc and clang). 
* Meson build system.
* Boost library.
* GNU Scientific Library.
* HDF5 library (Only C-bindings required, not the C++ API).
* Doxygen (only for documentation)
* Sphinx with Breathe extension (only for documentation)
* Python matplotlib (only for benchmark plots)

Building from Source
--------------------

The library is build using the modern Meson build system, which is 
available as standard package in most distributions and can also be 
installed via Python's pip. For more info, see
https://mesonbuild.com


To build the library::

  $> cd <repository>
  $> meson mbuild --buildtype=release
  $> cd mbuild
  $> ninja

This will compile with optimization and without debug symbols. Other
possibilities are `--buildtype=debug` and `--buildtype=debugoptimized`
To use a different compiler, e.g. clang, prefix the meson command
with `CC=clang CXX=clang++`.
See https://mesonbuild.com/Running-Meson.html for general Meson usage.


Installing
----------

To install the library systemwide::

  $> ninja install

To install in a user-defined location::

  $> DESTDIR=<prefix absolute path> ninja install


Using the Library
-----------------

The various header files needed for using the library are installed 
in a subdirectory `reprimand`. The executables have to be linked
with `libRePrimAnd`.


Creating Documentation
----------------------

To just build the documentation, use the target `documentation`::

  $> ninja documentation

The resulting pages can be found in the build directory under
`docsrc/sphinx/index.html`
When installing the library, the documentation is installed as well
by default, to `<prefix>/usr/local/share/doc/libreprimand/index.html`

The building of the documentation requires sphinx with the breathe 
extension as well as doxygen. To disable building documentation and 
remove the corresponding dependencies, use the build option::

  $> cd <repository>/mbuild
  $> meson configure -Dbuild_documentation=false


Running Tests
-------------

Please also take a minute to run the unit tests to ensure 
correct compilation::

  $> ninja test


Benchmarks
----------

The repository contains code to map the efficiency and accuracy of
the primitive recovery, producing the plots shown in the 
article. To recreate the data and plots::

  $> ninja benchplots
  $> ninja accuracyplots

The resulting pdf figures are placed in the build directory under
tests/benchmarks

This requires Python+matplotlib. To disable building benchmarking and 
remove the corresponding dependencies, use the build option::

  $> cd <repository>/mbuild
  $> meson configure -Dbuild_benchmarks=false

Visualizing Master Function
---------------------------

In addition, there is code to sample the primitive recovery master
function (the central ingredient of the scheme) for various cases,
as shown in the paper::

  $> ninja srootdata

The resulting data files are placed in the build directory under 
`tests/sample_root/`
 



