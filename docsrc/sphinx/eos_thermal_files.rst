EOS Files
=========

The framework contains functions to load thermal EOS from file.
The type of EOS is determined transparently to the user.
For this we defined a file format for thermal EOS, based on HDF5.

Loading EOS
-----------
Loading an EOS is simple, after including the header

.. code:: cpp

   #include "eos_thermal_file.h"
   
   using namespace EOS_Toolkit;

one can obtain an EOS via
   
.. code:: cpp

   units u = units::geom_solar();
   auto eos = load_eos_thermal("path/example.eos.h5", u);
  
  
The EOS files are based on SI units internally, therefore one
has to specify the geometric :cpp:class:`unit <EOS_Toolkit::units>`
system the returned EOS object 
should employ. In the example above, we use :math:`G=c=M_\odot=1`.

This version of the library comes with one example file
`EOS/HYB1.80_MS1_PP.eos.h5`, which is a hybrid EOS based on the MS1
cold EOS. 
 
 
Creating EOS Files
------------------
To create EOS files, the library provides a Python interface.
The EOS file format should not matter to users unless they want to
extend it with custom EOS types.
The interface consists of a single-file module
`EOS/create/reprimand_eos_format.py`.
The following EOS types are available:

Hybrid EOS
^^^^^^^^^^

.. py:currentmodule:: reprimand_eos_format

To create a hybrid EOS with cold part given by a tabulated barotropic EOS, use the following. See
:ref:`EOSColdTabulated` for details on the cold EOS.

.. autofunction:: save_thermal_hybrid_table

To create a hybrid EOS with cold part given by a piecewise polytropic EOS, use the following. See
:ref:`EOSPiecewisePoly` for details on the cold EOS.

.. autofunction:: save_thermal_hybrid_pwpoly

To create a hybrid EOS with cold part given by a polytropic EOS, use the following. See
:ref:`EOSPoly` for details on the cold EOS.

.. autofunction:: save_thermal_hybrid_poly

Ideal Gas EOS
^^^^^^^^^^^^^

.. autofunction:: save_thermal_idealgas
