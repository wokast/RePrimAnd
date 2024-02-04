.. _tov_solver_interface:

Interface
---------

There are two classes to represent spherical neutron stars. 
The most important is :cpp:class:`~EOS_Toolkit::spherical_star_properties`, 
which collects scalar measures and the EOS. A derived class,
:cpp:class:`~EOS_Toolkit::spherical_star`, additionally provides
the radial profiles describing metric and matter distribution.


There are different TOV solver functions for different use cases. 
One can either just extract global scalar measures, or also obtain the 
radial profiles found during the solution process. If profiles are not 
needed, one should use the first variant, which saves memory and 
computational costs since no interpolation tables have to be set up.

The function returning only global properties is called
:cpp:func:`~EOS_Toolkit::get_tov_properties` 
and the one returning everything is called
:cpp:func:`~EOS_Toolkit::get_tov_star`.

To control the accuracy, there is a class 
:cpp:class:`~EOS_Toolkit::star_accuracy_spec` which allows
to express the desired accuracies for mass, radius, 
moment of inertia, and tidal deformability individually, and
whether the tidal deformability and/or bulk radius is needed at all.
The specification of tolerances is completely independent from 
the solution method. The available TOV solver translates the
specification into suitable internal parameters, using heuristic
formulas that have been calibrated using the measured error for 
a large set of tabulated nuclear
physics EOS, polytropic EOS spanning a wide range of compactness, 
and piecewise polytropic EOS. 

There are two functions to create an error specification:
:cpp:func:`~EOS_Toolkit::star_acc_simple` and
:cpp:func:`~EOS_Toolkit::star_acc_detailed`.
The former distinguishes only the accuracy for deformability
from anything else, while the latter allows finer grained control.
The default values are appropriate for most applications. If very
high precision is needed one should narrow the tolerances, and
if speed is an issue one should widen the tolerance as much as 
possible. In particular, tidal deformability and bulk radius 
computation can be disabled if not needed in order to increase speed.


All solvers take EOS and central baryonic mass density as first arguments.
The third argument is the accuracy spec created by one of the above 
functions. 

The following example creates an EOS on the fly and computes a single TOV model

.. literalinclude:: minimal_tov.cc
   :language: cpp


The library also provides a method to find the maximum mass model
along a TOV sequence, :cpp:func:`~EOS_Toolkit::find_rhoc_tov_max_mass`.
Another function, :cpp:func:`~EOS_Toolkit::find_rhoc_tov_of_mass`,
finds a TOV model with given gravitational mass. Note this is not 
efficient for finding many models for a given EOS since it finds a
TOV solution for each step of the root finding. For such applications,
e.g. in GW data parameter estimation, it is much more efficient
to use the star sequence functionality provided by the library 
which is based on interpolation tables.


.. tip::

   NSs represented by  :cpp:class:`~EOS_Toolkit::spherical_star_properties` or
   :cpp:class:`~EOS_Toolkit::spherical_star` can be copied cheaply, and one 
   does not need to worry about memory management. 
   Internally, the EOS data and the NS profile are managed by reference counted
   shared pointers.


