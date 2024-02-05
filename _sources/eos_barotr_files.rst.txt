EOS Files
=========

The library defines a file format for barotropic EOS, based on HDF5.
Users should not need to deal with the format directly,
the framework contains functions to load and save barotropic EOS.
The type of EOS is determined while loading transparently to the user.

To avoid unit confusion, the files use SI units internally. When saving an EOS,
conversion happens automatically using the unit system stored within the EOS.
When loading, one has to specify the unit system, which is assumed to be geometric.

Note for developers of custom EOS:
The ability to load and/or save is optional for a given EOS type. Custom EOS only 
need to define those methods if needed. However, the library provides an abstracted 
interface for storing data without having to deal with the hdf5 interface directly.
This makes it easy to write file methods. As an example, see
`library/EOS_Barotropic/eos_barotr_spline_file.cc`.

Loading EOS
-----------
For loading a barotropic EOS from file, 
use :cpp:func:`~EOS_Toolkit::load_eos_barotr` as in the example below.

.. code:: cpp

   #include "eos_barotr_file.h"

   using namespace EOS_Toolkit;

   units u = units::geom_solar();
   auto eos = load_eos_barotr("path/example.eos.h5", u);
  
  
One has to specify the geometric :cpp:class:`unit <EOS_Toolkit::units>` 
system the returned EOS object should employ. In the example above, 
we use :math:`G=c=M_\odot=1`.

Example EOS
-----------

The library repository contains a small collection of EOS files 
representing barotropic nuclear physics EOS models. Most are 
spline-based (tabulated) and some are piecewise polytropic approximations.
The collection can be found in the folder
`EOS/nuclear_physics_eos_collection`. 

The original sources have been sanitized and resampled at our own 
discretion, and come with no guarantees. The collection is intended
for educational and testing purposes. 
 
Creating EOS Files
------------------

In order to create barotropic EOS to file, first create the EOS object and
the save it to a file using
:cpp:func:`~EOS_Toolkit::save_eos_barotr`.
For most use cases, it is probably more convenient to use the Python bindings, 
which include the C++ functions for creating and saving EOS. 



.. note::
   Older versions of the library used a separate Python module for creating EOS files,
   while the C++ interface could only load. This asymmetric design was abandoned,
   and the module removed in version 1.7. New EOS files should always be created
   using the regular Python interface to the C++ library.

.. note::
   Barotropic EOS of tabulated type are deprecated and
   cannot be saved. Existing EOS files for this type contained the original 
   irregularly spaced sample points from which such EOS had to be created by 
   resampling. Therefore, those EOS objects do not contain the information
   to save in the original format. The old format can still be read, however.
   It is recommended to use the new spline-based EOS instead.



