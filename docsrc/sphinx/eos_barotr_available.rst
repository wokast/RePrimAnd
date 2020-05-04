Available EOS
=============


Polytropic EOS
--------------

In terms of mass density :math:`\rho` or pseudo-enthalpy :math:`g`, 
the EOS is given by

.. math::
   :label: press_polytrope

   P &= \rho_p \left(\frac{\rho}{\rho_p}\right)^\Gamma 
     = \rho_p \left( \frac{g-1}{1+n} \right)^{1+n}  \\   
   \epsilon &= n \left(\frac{\rho}{\rho_p}\right)^\frac{1}{n} 
             = \frac{g-1}{\Gamma}
            
whith the adiabatic exponent :math:`\Gamma = 1 + \frac{1}{n}`.
The adiabtic index :math:`n>0` and the polytropic density 
scale :math:`\rho_p>0` are free parameters of the EOS.
The parameter :math:`\rho_p>0` is more commonly given in terms of 
the polytropic constant :math:`K=\rho_p^{1-\Gamma}`.
The drawback of the latter is that the unit depends on :math:`\Gamma`.

It is easy to show that the EOS is isentropic, and we further define it
as a zero-temperature EOS.
Zero-temperature polytropes are 
used frequently as a toy model for degenerate neutron star matter
and for testing numerical relativity codes. 

The pseudo-enthalpy is simply :math:`g=h/h_0`. Further, the minimum
enthalpy is :math:`h_0=1`.
The following relations apply:

.. math::
   g(\rho) &=  1 + (n+1) \left(\frac{\rho}{\rho_p}\right)^{1/n} \\
   \rho(g) &=  \rho_p \left( \frac{g-1}{1+n} \right)^n \\
   c_s^2   &=  \frac{g-1}{ng} 
   
This EOS does not provide an electron fraction (in contrast to more
realistic on nuclear physics models which compute beta-equilibrium).

.. warning::

   If :math:`n<1` the soundspeed would exceed the speed of light 
   above a critical density. In this case, the validity range specified 
   by the user is therefore automatically adjusted to ensure

   .. math::

      g < \frac{1}{1-n}

Piecewise Polytropic EOS
------------------------

The piecewise polytropic EOS consists of polytropic segments with
different polytropic exponents and an additional offset of specific 
energy. 

.. math::

   P(\rho) 
    &= \rho_{p,i} \left(\frac{\rho}{\rho_{p,i}}\right)^{\Gamma_i}
    \qquad \mathrm{for}\quad \rho_{b,i} \le \rho < \rho_{b,i+1} \\
   \epsilon(\rho) 
     &= n_i \left(\frac{\rho}{\rho_{p,i}}\right)^\frac{1}{n_i} 
        + \epsilon_i

where :math:`\rho_{b,i}` denote the mass density at the segment 
boundaries, :math:`\rho_{p,i}` and :math:`\Gamma_i=1+1/n_i` 
the polytropic 
density scale and polytropic exponent of segment :math:`i`.
Compared to regular polytropes, each segment has an additional
constant offset :math:`\epsilon_i` for the specific energy.

The EOS is completely determined by the exponents :math:`\Gamma_i`,
the segment boundaries :math:`\rho_{b,i}`, and :math:`\rho_{p,0}`
of the first segment. Matching pressure and energy across boundaries 
then determines :math:`\epsilon_{i>0}` and :math:`\rho_{p,i>0}`.
Further, we set :math:`\epsilon_0=0`.

It is easy to show that the EOS is isentropic, and we further define it
as a zero-temperature EOS.

The minimum enthalpy is :math:`h_0=1`.
Pseudo-entaply and enthalpy are identical, i.e. :math:`g=h`. 
In terms of the pseudo-enthalpy, the segments are given by

.. math::
   g(\rho) 
     &= 1 +  \epsilon_i
          + (n_i+1) \left(\frac{\rho}{\rho_{p,i}}\right)^{1/n_i} \\
   \rho(g) 
     &= \rho_{p,i} \left( \frac{g-1-\epsilon_i}{1+n_i} \right)^{n_i} \\
   P(g) 
    &= \rho_{p,i} \left( \frac{g-1-\epsilon_i}{1+n_i} \right)^{1+n_i} \\
   \epsilon(g)
    &= \frac{g-1-\epsilon_i}{\Gamma_i} + \epsilon_i \\
   c_s^2(g) 
    &= \frac{g-1-\epsilon_i}{n_ig}

.. warning::
   Segments with :math:`n_i<1` would lead to superluminal soundspeed,
   if the corresponding critical density falls within the segment.
   The user-specified validity range is automatically reduced to 
   prevent this, if necessary.

.. warning::
   The intended use case for this EOS contains just few (<10) segments.
   It would be very inefficient to approximate arbitrary EOS using 
   hundreds of segments. For this, use the tabulated EOS below.

Tabulated EOS
-------------

This is the most general EOS, where all functions are implemented
as efficient linear lookup tables. Those lookup tables are regularly 
spaced, hence evaluation cost is independent of table resolution.
The low convergence order of linear lookup can be compensated by 
table resolution.

To create a tabulated EOS, one provides vectors consisting of data 
at arbitrarily spaced sample points

The provided sample points are interpolated
to regularly spaced samples for the lookup tables using  
monotonic cubic spline interpolation. This avoids violation of
monotonicity conditions by overshoots. 

Since the lookup tables have to cover orders of magnitude, logarithmic 
spacing is employed for the independent variable. However, a constant 
offset is added to prevent wasting many sample points on very low 
values, which could otherwise happen if the lowest sample is very close 
to zero. 
The heuristic algorithm determines the dynamic range spanned from the 
median to the maximum value of the mass density, and 
adjusts the offset such that this range is covered by half of the 
lookup table points. The lookup tables in terms of pseudo-enthalpy 
:math:`g-1` adjust the offset such that the shifted value spans the 
same magnitude range as for the density.

Most functions are tabulated with the pseudo enthalpy :math:`g` as
independent variable, in addition the mapping between :math:`g` and
:math:`\rho` is tabulated.

Temperature and electron fraction can be provided optionally.
If the temperature is not provided when creating the EOS, it is assumed
to be a zero-temperature EOS.

Finally, between zero density and the lowest density covered by the 
lookup table, a polytropic EOS with an additional offset in specific
energy is attached. The polytropic exponent is a free parameter, offset
and polytropic density scale are fixed my matching conditions.
Temperature and electron fraction, if available, are kept constant.
The motivation is that nuclear physics EOS tables typically do not 
extend down to zero density, but for setting up intial data, such as
neutron stars, it is convenient not to have to worry about artificial
cutoffs before reaching the surface. For use in evolution codes, it 
is up to the user to make sure results do not rely on this extension,
e.g. highly diluted expelled matter. 

.. warning::
   Due to the heuristic setup of lookup tables and the required 
   interpolation, the provided sample points should have a roughly 
   uniform coverage of logarithmic density between median value and 
   maximum.













