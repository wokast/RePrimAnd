EOS Files
=========

The framework contains functions to load thermal EOS from file.
The type of EOS is determined transparently to the user.
For this we defined a file format for thermal EOS, based on HDF5.

Loading EOS
-----------
For loading a thermal EOS from file, use 
:cpp:func:`~EOS_Toolkit::load_eos_thermal` as in the example below.

.. code:: cpp

   #include "eos_thermal_file.h"
   
   using namespace EOS_Toolkit;

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



In order to create barotropic EOS to file, first create the EOS object and
the save it to a file using
:cpp:func:`~EOS_Toolkit::save_eos_thermal`.
For most use cases, it is probably more convenient to use the Python bindings, 
which include the C++ functions for creating and saving EOS. 




.. note::
   Older versions of the library used a separate Python module for creating EOS files,
   while the C++ interface could only load. This asymmetric design was abandoned,
   and the module removed in version 1.7. New EOS files should always be created
   using the regular Python interface to the C++ library.
