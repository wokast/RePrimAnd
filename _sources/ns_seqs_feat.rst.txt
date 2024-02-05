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

The framework provides a representation of generic sequences and
a specialization representing stable branches.
Sequences are spanning a range of central densities, and present
NS properties as functions of central pseudo-enthalpy.
Branches are sequences that restrict the range to a 
segment where the gravitational mass is strictly increasing with 
pseudo-enthalpy, and which extends up to the maximum mass. 

Branches also offer all properties as function of the gravitational 
mass, and they provide the maximum mass as wel as the properties of
the corresponding NS model.


We note that the EOS framework insists on causality, and 
therefore the same guarantee extends to the star sequences. By design,
the library actively prevents the computation of TOV solutions with 
causality-violating regions. Even though the equations for the static 
TOV solution do not break down for such EOS, this does not hold
for dynamic evolution equations that govern, e.g., oscillations. 
Such models are considered physically meaningless.


It may happen that the central density exceeds the EOS validity
range before the maximum NS mass is reached. The library
distinguishes this case, and allows to query whether a branch
extends to the physical maximum or the one given by the EOS range.

