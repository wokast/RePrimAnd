Features
--------

The library provides a method to recover primitive variables 
:math:`\rho, \epsilon, Y_e, P, v^i, W, E^i` from the evolved 
variables :math:`D, \tau, S_i, Y_e^T, B^i` of the conservation-law
formulation of relativistic ideal magnetohydrodynamics.
The code implements exactly the algorithm described in the acompanying
article.

The implementation determines whether the evolved variables 
correspond to valid primitives or not. Valid means that the matter 
state :math:`\rho,\epsilon,Y_e` is in the valid range of the EOS, in 
particular that the fluid internal energy is above the value at zero 
temperature. In addition, one can specify technical restrictions if 
meaningful for the evolution code: a speed limit and bound for 
magnetization.

Further, the code can project the evolved variables back onto
the valid regime. The main correction is to set the fluid internal 
energy to the closest value valid at given mass density for the EOS.
In addition, the velocity can be limited to the speed limit by 
rescaling the momentum, and the electron fraction can be adjusted to 
the valid range of the EOS. The magnetic field is never corrected,
as this cannot be done locally in a meaningful way.

If a given correction is applied is decided by a parametrized error
policy. The implemented policy is geared towards simulations of
isolated neutron stars, binary neutron mergers, and collapse to black 
hole. It only contains simple conditions, as everything more 
involved is likely specific to the evolution code.

In detail, the policy is characterized by the following parameters:

1. A density above which the only allowed corrections are adjusting 
   internal energy below zero temperature value and the electron 
   fraction to the allowed range ("strict" regime).
2. A maximum velocity above which failure is reported in the strict 
   regime and momentum rescaling applied otherwise. 
3. An upper limit for the ratio of magnetic and mass energy densities,
   above which failure is reported (at any density).
4. A flag whether to restrict the electron fraction to the 
   allowed range also in the strict regime, or whether to report 
   failure.

It is simple to declare different error policies and apply them 
in different regimes. For example, the evolution code could use the 
lapse function to estimate if a point is inside a black hole, and 
then apply a more lenient policy.

Finally, the implementation can enforce an artificial atmosphere on the
fluid. In detail, below a given mass density, :math:`\rho,\epsilon,Y_e`
are set to specified values, and the velocity set to zero. For the EM 
part, the electric field is set to zero and the magnetic field is left 
unchanged. The evolved variables are set consistent with those 
primitives. This is a very crude approach that might lead to incorrect 
evolution of the magnetic field in regions with atmosphere. Any
improvement is specific to the evolution code and should be handled 
there.


