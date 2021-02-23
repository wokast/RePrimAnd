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
* Python pybind11 package (only for Python bindings)

Building from Source
--------------------

The library is build using the modern
`Meson <https://mesonbuild.com>`_
build system, which is 
`available <https://mesonbuild.com/Getting-meson.html>`_
as standard package in most distributions and can also be installed 
via Python's pip.

To build the library, 

.. code::

   cd <repository>
   meson mbuild --buildtype=release --prefix=<custom install location>
   cd mbuild
   ninja
   
This will compile with optimization and without debug symbols. Other
possibilities are `--buildtype=debug` and `--buildtype=debugoptimized`
To use a different compiler, e.g. clang, prefix the meson command
with `CC=clang CXX=clang++`.
See `here <https://mesonbuild.com/Running-Meson.html>`_ for general 
Meson usage.


Installing
^^^^^^^^^^

To install the library, use

.. code::

   cd mbuild
   ninja install
   
This will install in a user-defined location if the `--prefix` option
was given during the build setup, otherwise systemwide. 


Using the Library
^^^^^^^^^^^^^^^^^

The various header files needed for using the library are installed 
in a subdirectory `reprimand`. The executables have to be linked
with library `RePrimAnd`.

A minimal example can be found in `example/minimal.cc`. Assuming
the library is installed where the compiler (below we use gcc as 
example) can find it, compilation should be straightforward:

.. code:: bash

   g++ -lRePrimAnd --std=c++11 minimal.cc


Python Bindings
---------------

The Python bindings are automatically build together with the library
itself, provided that a Python3 installation is detected which contains
the pybind11 package. The latter can be easily installed via pip or 
conda. To disable building the Python bindings and remove the 
corresponding dependencies, use the build option

.. code::

   meson configure -Dbuild_python_api=false


If the Python bindings have been build, they are automatically installed 
together with the library. The Python extension module is called 
`pyreprimand`.


Creating Documentation
----------------------

To just build the documentation, use the target `documentation`.

.. code::
 
   ninja documentation

The resulting pages can be found in the build directory under
`docsrc/sphinx/index.html`.

The building of the documentation requires sphinx with the breathe 
extension as well as doxygen (Note sphinx breathe currently 
has problems with the C++11 trailing
return type syntax, misreporting return types as `auto`.)

To disable building documentation and remove 
the corresponding dependencies, use the build option

.. code::

   meson configure -Dbuild_documentation=false


If the documentation is build, it is installed automatically when 
installing the library, by default to 
`<prefix>/usr/local/share/doc/libreprimand/index.html`.


Running Tests
-------------

Unit tests
^^^^^^^^^^

Please also take a minute to run the unit tests to ensure 
correct compilation

.. code::

   ninja test
   

Benchmarks
^^^^^^^^^^

The repository contains code to map the efficiency and accuracy of
the primitive recovery, producing the plots shown in the 
article. To recreate the data and plots,


.. code::

   ninja benchplots
   ninja accuracyplots
   
The resulting pdf figures are placed in the build directory under
`tests/benchmarks`.

This requires Python+matplotlib. To disable building benchmarking and 
remove the corresponding dependencies, use the build option

.. code::

   cd <repository>/mbuild
   meson configure -Dbuild_benchmarks=false

Visualizing Master Function
^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition, there is code to sample the primitive recovery master
function (the central ingredient of the scheme) for various cases,
as shown in the paper.

.. code::

   ninja srootdata
  
The resulting data files are placed in the build directory under 
`tests/sample_root/`.
 



