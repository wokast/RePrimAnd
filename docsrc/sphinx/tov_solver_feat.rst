Features
========

The library contains functions for computing properties of spherical neutron stars.
It solves the TOV ODE as well as the ODE for tidal deformability. Models are specified
by central density and the EOS (via the barotropic EOS interface).
Basic support for finding models by mass or the maximum mass model
is also implemented.

Any barotropic EOS can be used, including ones which are not adiabatic (for example,
density-dependent arbitrary electron fraction or non-zero temperature). 
One can also employ EOS with phase transitions where the pressure remains constant 
over a density range and the soundspeed drops to zero.
The TOV solver employs a formulation of the TOV equations suitable also for those cases.

The following global properties are computed:

* Gravitational (ADM) mass :math:`M_g`
* Baryonic mass :math:`M_b`
* Circumferential surface radius :math:`R_c`
* Proper volume :math:`V_p`
* Moment of inertia :math:`I` 
* Central internal energy :math:`\epsilon_c`, pressure :math:`P_c`, and 
  sound speed :math:`c_{sc}` 
* Dimensionless tidal deformability :math:`\Lambda` and love number :math:`k_2`
* Properties of the "bulk" (optional, see below) 

The metric can also be obtained. It is expressed as follows

.. math::
   \mathrm{d} s^2 = - e^{2\nu(r_c)} \mathrm{d} t^2 + e^{2\lambda(r_c)} \mathrm{d} r^2 
                     + r_c^2 \mathrm{d} \Omega^2 

The metric potentials :math:`\nu(r_c)` and :math:`\lambda(r_c)` 
and the matter state are provided as 
functions of circumferential radius :math:`r_c`.

The tidal deformability formalism is described in 
:footcite:p:`Hinderer:2008:1216, Han:2019:083014`.
Since the formalism assumes adiabatic perturbations, the tidal deformability is 
only computed for isentropic EOS, for example zero-temperature stars. 
The library uses a reformulation of the ODE in :footcite:p:`Han:2019:083014` that is well-behaved across phase transitions (currently undocumented).

The moment of inertia is defined in the slow-rotation approximation following 
:footcite:p:`Hartle:1967:1005`. 
Finally, the library can also compute a relatively new definition of "bulk" properties 
from :footcite:p:`Kastaun:2016`. This definition uses a sort of maximum compactness 
iso-density surface and can also be applied to systems without symmetry or surface, 
such as early merger remnants.

.. warning::

    The TOV solver functionality in general is experimental.
    In particular, phase transitions are allowed by design but they 
    are currently not tested at all.

.. note::

   Since tabulated EOS and piecewise polytropic EOS are not 
   differentiable at segment boundaries, for those EOS the convergence 
   speed of ODE
   solvers effectively drops to first order with increasing resolution. 
   If speed is an issue it is recommended not to request solutions 
   with relative errors below :math:`10^{-10}`. 
   A fully differentiable tabulated EOS is planned.

   
References
^^^^^^^^^^

.. footbibliography::

