.. _python_interface:

Python Bindings
---------------

The `RePrimAnd` library also provides a Python3 extension module 
`pyreprimand` that allows the use of the library from within
Python.


As an added benefit, functions from the C++ interface are vectorized
in the Python interface whenever it makes sense, meaning that it is 
possible to evaluate those functions for `numpy` 
arrays of arbitrary dimension. This should be very efficient since the 
loops are performed in compiled C++ code.

EOS Objects
^^^^^^^^^^^

The Python bindings include the functionality of barotropic and thermal
EOS, as well as the functions to load and save these. 
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
   
   u   = pyr.units.geom_solar()
   eos = pyr.load_eos_barotr("EOS/MS1_PP.eos.h5", units=u)
   
   rho = np.linspace(eos.range_rho.min, eos.range_rho.max, 1000)
   press = eos.press_at_rho(rho)
   
Creating EOS Files
^^^^^^^^^^^^^^^^^^

One prime use for the Python interface is the creation of EOS files based on
data from other sources. Polytropic, piecewise polytropic, and spline-based 
EOS for tabulated data are supported directly. Other EOS types can be represented 
by first sampling them.

Creating a piecewise polytropic EOS file is as simple as

.. code-block:: python

   u   = pyr.units.geom_solar()

   eos = pyr.make_eos_barotr_pwpoly(rho_poly, 
                 rho_bounds, gammas, rho_max, uc)

   pyr.save_eos_barotr("example_pp.eos.h5", eos)

The arguments `rho_bounds` and `gammas`, which would be vectors in the C++ interface,
can be passed as numpy arrays in Python. The dimensionful arguments
`rho_poly` and `rho_bounds` have to be specified in the unit system desired 
for the EOS, specified by the `uc` argument with respect to SI units.
See also :cpp:func:`~EOS_Toolkit::make_eos_barotr_pwpoly`.

For creating an EOS file that represents tabulated data, use the spline-based
EOS as follows

.. code-block:: python

   u   = pyr.units.geom_solar()

   eos = pyr.make_eos_barotr_spline(rho, eps, press, csnd, 
                temp, efrac, is_isentropic, range_rho, n_poly, 
                eos_units, pts_per_mag)

   pyr.save_eos_barotr("example_spline.eos.h5", eos)

where `rho`, `eps`, `press`, `csnd`, `temp`, and `efrac` are numpy 
arrays with the arbitrarily-spaced tabulated sample points. If temperature 
and/or electron fraction are not available, pass an empty list instead.
Any dimensionful arguments have to be specified in the unit system desired 
for the EOS, specified by the `uc` argument with respect to SI units.
See also :cpp:func:`~EOS_Toolkit::make_eos_barotr_spline` 


.. note:
   There used to be a seperate Python module for creating EOS files. This is deprecated now.
   The old tabulated EOS type is also deprecated and cannot be saved via the current
   library. Use spline-based EOS above instead.




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
will not work. See below for working with TOV sequences.


NS Sequences
^^^^^^^^^^^^

The interface for NS sequences and their (stable) branches is available 
in Python. In particular, use the function 
:cpp:func:`~EOS_Toolkit::make_tov_branch_stable` for computing a stable branch
for TOV solutions, the function :cpp:func:`~EOS_Toolkit::save_star_branch`
to save it to file, and the function :cpp:func:`~EOS_Toolkit::load_star_branch`
to load such a file.

The methods of objects :cpp:class:`~EOS_Toolkit::star_branch` (representing stable 
branches) are vectorized, i.e., they accept `numpy` arrays. If the input is outside 
the valid range, the result is NAN. The same holds for :cpp:class:`~EOS_Toolkit::star_seq` 
objects representing a general NS sequence.

An example Python script plotting TOV sequences can be found under
`examples/pwpoly_TOV.py`.

Units
^^^^^

The Python bindings include the :cpp:class:`EOS_Toolkit::units` 
objects, which are used throughout the interface to specify unit systems
in a unified consistent way. The unit objects should also be very
useful on their own, for example in quick interactive calculations.
The very popular system defined by :math:`G=c=M_\odot=1` is predefined as
`pyr.units.geom_solar()`. Note this is a function; when setting up
geometric units, the default values for :math:`G` and/or :math:`M_\odot` 
can be overridden by parameters msun_si, g_si. 
Unit objects provide specific units via data members such as 
`u.density` (try Tab-completion for the complete list!).
Note that temperatures in the EOS framework are always in `MeV` and the unit
conversion objects do not define a temperature unit. 
