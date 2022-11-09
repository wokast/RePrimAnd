# What is this?

The pyRePrimAnd package provides a Python interface
for the RePrimAnd C++ library which can be found 
[here](https://wokast.github.io/RePrimAnd/index.html).

RePrimAnd is a support library for numerical simulations of general 
relativistic magnetohydrodynamics and other neutron star related
problems. RePrimAnd provides 

* A general framework for handling matter equations of state, 
* A solver for the TOV equations describing nonrotating neutron stars
* Tools for precomputing properties for sequences of neutron stars.
* Methods for recovering primitive variables in GRMHD. This is not 
  included in the Python interface because it is mainly needed for 
  high performance computing.

## Documentation

The documentation for the library and the Python interface 
can be found [here](https://wokast.github.io/RePrimAnd/index.html)

## Installation

For the Linux platform, we provide Python wheels bundled with the 
precompiled library. They can be installed from 
[pypi](https://pypi.org) using


```bash
pip install pyreprimand
```

For other platforms, including Macs, one would first have to build
and install the reprimand library from 
source (see [here](https://github.com/wokast/RePrimAnd)).
Installing pyreprimand with pip will then try building from
source. This may still fail as MacOS is not the development platform.
In this case, please open a ticket on the issue tracker


## Requirements

* The numpy package

Only when building from source distribution:

* A C++11 capable compiler (tested with gcc and clang). 
* Python pybind11 package >= 2.6.0 (only for Python bindings)
* The RePrimAnd library (only when building from source dist)

## Support

In case of errors in the Python interface or the library, submit a ticket to
the [issue tracker](https://github.com/wokast/RePrimAnd/issues).



