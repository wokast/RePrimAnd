.. _python_interface:

Python Bindings
---------------

The `RePrimAnd` library also provides a Python3 extension module 
`pyreprimand` that allows the use of the EOS framework from within
Python3. The main purpose of this is postprocessing simulation data. 
It is therefore possible to evaluate the EOS methods for `numpy` 
arrays of arbitrary dimension. This should be very efficient since the 
loops are performed in compiled C++ code.

The Python bindings include the functionality of barotropic and thermal
EOS, as well as the functions to load these from file. It does not 
expose the low-level infrastructure, such as EOS implementation objects,
or direct constructors. In addition, the Python bindings include the 
:cpp:class:`unit` objects, which are needed when loading an EOS, but 
should be useful in general for unit conversion.

The Python bindings expose nearly all EOS methods, using the same naming
as the C++ interface. One difference is that the Python bindings do not 
include the methods based on matter state objects, as this interface 
is not suitable for work with arrays. Only the C++ methods for direct
evaluation of a single quantity, such as 
:cpp:func:`~eos_thermal::press_at_rho_eps_ye` 
are offered. As with the C++ pendants, they return NAN for invalid 
input. Range-checking methods such as 
:cpp:func:`~eos_thermal::is_rho_valid()` are also 
exposed and can be used with numpy arrays, resulting in `numpy` arrays
of type `bool`. Another minor difference is that C++ getter methods 
such as :cpp:func:`~eos_thermal::range_rho` are mapped to Python 
object properties, so in Python one should not add brackets ().

The example code below shows how to load a barotropic EOS and evaluate 
the pressure for a numpy array of mass densities spanning the full EOS
valid range.

.. code-block:: python

   import pyreprimand as pyr
   import numpy as np
   
   u   = pyr.units.geom_solar
   eos = pyr.load_eos_barotr("EOS/MS1_PP.eos.h5", units=u)
   
   rho = np.linspace(eos.range_rho.min, eos.range_rho.max, 1000)
   press = eos.press_from_rho(rho)
   


