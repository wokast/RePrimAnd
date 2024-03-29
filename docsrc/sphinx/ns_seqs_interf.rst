.. _ns_seqs_interface:

Interface
---------

There are two classes to represent neutron star sequences. 
The basic one is :cpp:class:`~EOS_Toolkit::star_seq`, 
which represents sequences with an arbitrary range of central density.
The other is :cpp:class:`~EOS_Toolkit::star_branch`, which represents 
a stable branch where properties can be obtained also as function of 
gravitational mass. The class :cpp:class:`~EOS_Toolkit::star_branch`
derives from :cpp:class:`~EOS_Toolkit::star_seq` and can thus be
used as a sequence as well. Both classes store the units they use,
which can be freely specified during creation but are assumed to be 
geometric. Each object provides a method
:cpp:func:`~EOS_Toolkit::star_seq::units_to_SI` to query the units.

To compute a sequence of TOV stars, use the function
:cpp:func:`~EOS_Toolkit::make_tov_seq`, specifying
the EOS, the desired accuracy (see TOV solver), 
the range of central pseudo-enthalpy, and the number 
of sample points. The units are taken from the EOS.

To compute a stable branch for TOV stars, use the function
:cpp:func:`~EOS_Toolkit::make_tov_branch_stable`.
To select the correct branch in case there are several, one can provide 
one central pseudo-enthalpy within the desired branch. 
The default value should work for almost any remotely realistic NS EOS. 
One can also specify a minimum mass that has to be covered by 
the sequence. The units for the sequence are copied from the EOS.

To create a star sequence directly from data points,
use :cpp:func:`~EOS_Toolkit::make_star_seq`, specifying 
vectors for the NS properties. Currently, those need to
be equally spaced in central pseudo enthalpy. For the latter,
only the range needs to be given, not a vector. Further,
one needs to specify the unit system to which the provided 
data refers, which is assumed geometric.


The following example demonstrates basic use of stable TOV branches.


.. literalinclude:: minimal_seq.cc
   :language: cpp



.. tip::

   NS sequences represented by  :cpp:class:`~EOS_Toolkit::star_sequence` or
   :cpp:class:`~EOS_Toolkit::star_branch` can be copied cheaply, and one 
   does not need to worry about memory management. 
   Internally, the data points for interpolation are managed by reference counted
   shared pointers. Using sequences or branches is thread-save.


