EOS Files
=========

The framework contains functions to load barotropic EOS from file.
The type of EOS is determined transparently to the user.
For this we defined a file format for barotropic EOS, based on HDF5.

Loading EOS
-----------
Loading an EOS is simple, after including the header

.. code:: cpp

   #include "eos_barotr_file.h"

   using namespace EOS_Toolkit;

one can obtain an EOS via
   
.. code:: cpp

   units u = units::geom_solar();
   auto eos = load_eos_barotr("path/example.eos.h5", u);
  
  
The EOS files are based on SI units internally, therefore one
has to specify the geometric :cpp:class:`unit <EOS_Toolkit::units>` 
system the returned EOS object 
should employ. In the example above, we use :math:`G=c=M_\odot=1`.

This version of the library comes with one example file
`EOS/MS1_PP.eos.h5`, which is a tabulated EOS created from a piecewise 
polytropic approximation of the `MS1` EOS. 

 
Creating EOS Files
------------------
To create EOS files, the library provides a Python interface.
The EOS file format should not matter to users unless they want to
extend it with custom EOS types.
The interface consists of a single-file module 
`EOS/create/reprimand_eos_format.py`.
The following EOS types are available:


.. py:currentmodule:: reprimand_eos_format

To create a tabulated barotropic EOS file, use the following. See
:ref:`EOSColdTabulated` for details on the EOS.

.. autofunction:: save_barotr_poly


To create a piecewise polytropic EOS file, use the following. See
:ref:`EOSPiecewisePoly` for details on the EOS.

.. autofunction:: save_barotr_pwpoly

To create a polytropic EOS file, use the following. See
:ref:`EOSPoly` for details on the EOS.

.. autofunction:: save_barotr_table
