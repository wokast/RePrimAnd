Features
========

The library contains functions for precomputing properties for sequences of 
neutron stars and interpolating those to any parameter. Further, the sequences
can be written and loaded from files in a hdf5 based format. Sequences can be 
computed using the TOV solver provided by the library, but also created using 
sampled data from other sources (note this would also allow for NS solutions from
modified GR theories, provided that the same properties are still meaningful).
The stored NS properties are: gravitational mass, baryonic mass, central matter 
state, moment of inertia, circumferential radius, and tidal deformability. 

The framework distinguishes generic sequences specified by central matter state,
and stable branches which are segments where the gravitational mass is strictly
increasing. Sequences are parametrized by central pseudo-enthalpy. Branches are
also parametrized as function of gravitational mass. They also provide the maximum 
mass model properties. Finally, they distinguish the corner case where the maximum 
mass model is given by the EOS validity range from a physical maximum of the TOV
solutions. Since the EOS framework insists on causality, the same guarantee
extends to the star sequences.

.. warning::

    The NS sequence functionality in general is experimental.




