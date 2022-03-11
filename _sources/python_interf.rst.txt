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
:cpp:class:`EOS_Toolkit::units` objects, which are needed when loading an EOS, but 
should be useful in general for unit conversion.

The Python bindings expose nearly all EOS methods, using the same naming
as the C++ interface. One difference is that the Python bindings do not 
include the methods based on matter state objects, as this interface 
is not suitable for work with arrays. Only the C++ methods for direct
evaluation of a single quantity, such as 
:cpp:func:`~EOS_Toolkit::eos_thermal::press_at_rho_eps_ye` 
are offered. As with the C++ pendants, they return NAN for invalid 
input. Range-checking methods such as 
:cpp:func:`~EOS_Toolkit::eos_thermal::is_rho_valid()` are also 
exposed and can be used with numpy arrays, resulting in `numpy` arrays
of type `bool`. Another minor difference is that C++ getter methods 
such as :cpp:func:`~EOS_Toolkit::eos_thermal::range_rho` are mapped to Python 
object attributes, so in Python one should not add brackets ().

The example code below shows how to load a barotropic EOS and evaluate 
the pressure for a numpy array of mass densities spanning the full EOS
valid range.

.. code-block:: python

   import pyreprimand as pyr
   import numpy as np
   
   u   = pyr.units.geom_solar
   eos = pyr.load_eos_barotr("EOS/MS1_PP.eos.h5", units=u)
   
   rho = np.linspace(eos.range_rho.min, eos.range_rho.max, 1000)
   press = eos.press_at_rho(rho)
   


TOV Solver
^^^^^^^^^^

The C++ TOV solver functionality is also accessible from Python. The functions 
:cpp:func:`~EOS_Toolkit::make_tov_star` and :cpp:func:`~EOS_Toolkit::get_tov_star_properties`
allow the creation of a single star with and without profile information, respectively. 
Constructing the objects :cpp:struct:`~EOS_Toolkit::tov_acc_simple` or 
:cpp:struct:`~EOS_Toolkit::tov_acc_precise` describing the desired accuracy 
works in Python as in C++.

A minor difference to C++ is that getter methods 
such as :cpp:func:`~EOS_Toolkit::spherical_star_properties::grav_mass` are mapped to Python 
object attributes, so in Python one should not add brackets ().
Also, since lambda is a reserved word in Python, the tidal deformability
is called `tov.deformability.lambda_tidal` compared to
`tov.deformability().lambda` in C++.

Methods related to stellar profiles, such as 
:cpp:func:`~EOS_Toolkit::spherical_star::press_from_rc`  
are vectorized and accept a numpy array as argument for the radius.
Note that computing TOV solutions is not vectorized, so passing a numpy array as central density
will not work. Future versions may contain dedicated functions for TOV sequences.

An example Python script plotting TOV sequences can be found under
`examples/pwpoly_TOV.py`.


 

