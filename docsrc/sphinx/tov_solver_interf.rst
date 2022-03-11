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
:cpp:func:`~EOS_Toolkit::get_tov_star_properties` 
and the one returning everything is called
:cpp:func:`~EOS_Toolkit::make_tov_star`.

Furthermore, there are two options regarding the accuracy. 
The fast option uses adaptive step size control, trying to archive 
given tolerances but does not guarantee meeting them.
A second option repeatedly employs the first one with increasingly 
small tolerances until changes fall below given accuracy.

All solvers take EOS and central baryonic mass density as first arguments.
The third argument is the accuracy. Passing an object of type 
:cpp:struct:`~EOS_Toolkit::tov_acc_simple` allows to specify tolerances
for the fast adaptive solution. Passing an object of type 
:cpp:struct:`~EOS_Toolkit::tov_acc_precise` instead allows to prescribe absolute 
error bounds for individual properties, such as mass, radius, and tidal deformability.

Further optional arguments specify whether to compute the "bulk" quantities and the 
tidal deformability. If not needed, computation of those can be turned off to increase 
performance.

The following example creates an EOS on the fly and computes a single TOV model

.. literalinclude:: minimal_tov.cc
   :language: cpp


The library also provides a method to find the maximum mass model
along a TOV sequence, :cpp:func:`~EOS_Toolkit::find_rhoc_tov_max_mass`.
Another function, :cpp:func:`~EOS_Toolkit::find_rhoc_tov_of_mass`,
finds a TOV model with given gravitational mass. Note this is not 
efficient for finding many models for a given EOS since it finds a
TOV solution for each step of the root finding. For such applications,
e.g. in GW data parameter estimation, using an interpolation table 
for a TOV sequence would be much faster.


.. tip::

   NSs represented by  :cpp:class:`~EOS_Toolkit::spherical_star_properties` or
   :cpp:class:`~EOS_Toolkit::spherical_star` can be copied cheaply, and one 
   does not need to worry about memory management. 
   Internally, the EOS data and the NS profile are managed by reference counted
   shared pointers.


