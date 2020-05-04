Available EOS
=============

Classical Ideal Gas
-------------------
This EOS serves mainly as a toy model useful for debugging.
It is given by

.. math::
   :label: press_idealgas

   P (\rho,\epsilon,Y_e) = (\Gamma - 1) \rho \epsilon
   
The adiabatic exponent :math:`\Gamma = 1 + 1/n` with :math:`n>0` 
is a free parameter. The pressure is independent of 
electron fraction, which is a dummy variable for this EOS.  
The validity range is given by

.. math::

   0 \le \rho \le \rho_\mathrm{max} \\
   0 \le \epsilon \le \epsilon_\mathrm{max} \\
   0 \le Y_e \le 1

where :math:`\rho_\mathrm{max}` and :math:`\epsilon_\mathrm{max}`
are free parameters specified by the user.

The adiabatic soundspeed is given by

.. math::

   c_s^2 = \frac{\left(\Gamma - 1\right) \epsilon}{\epsilon + 1/\Gamma} 

In case that :math:`\Gamma > 2`, the following restriction is
automatically applied in order to ensure subluminal speed of sound 
and dominant energy condition. 

.. math:: 

   \epsilon_\mathrm{max} \le \frac{1}{\Gamma (\Gamma - 2)}



For the classical ideal gas law, :math:`T \propto \epsilon`, and zero 
temperature implies zero pressure.
However, this law breaks down near zero temperature since it 
leads to a divergence of the entropy. One can avoid this by using a
different relation for the temperature, which is however not unique.
Consequently, our implementation
of this EOS does not provide temperature and specific entropy.

.. note::

   Since the pressure is zero at zero temperature, a classical ideal gas 
   is not useful for neutron stars and such. One could keep 
   :math:`P (\rho,\epsilon,Y_e)` but redefine temperature to avoid 
   zero pressure, emulating a degeneracy pressure. Since the adiabatic 
   curves for Eq. :math:numref:`press_idealgas` are given by polytropes, 
   :math:`P` and :math:`\epsilon` at zero temperature would follow a 
   polytrope as well. The resulting EOS is what is usually called ideal 
   gas in numerical relativity. To get such an EOS, use the 
   :ref:`def_hybrid_eos` with a polytrope as cold EOS. 
   The difference is simply the lower bound for :math:`\epsilon`,
   which becomes a strictly positive function of density.
   
   
.. _def_hybrid_eos:

Hybrid EOS
----------
This EOS is based on an barotropic EOS for zero temperature
together with a simple analytic prescription for thermal contributions.

.. math::
   :label: def_hybrideos
   
   P (\rho,\epsilon,Y_e) 
   = P_c(\rho) 
     + (\Gamma_\mathrm{th} - 1) \rho (\epsilon - \epsilon_c(\rho))

where :math:`P_c(\rho)` and :math:`\epsilon_c(\rho)` defines the 
zero-temperature EOS.
The pressure is independent of 
electron fraction, which is a dummy variable for this EOS.  
The validity range is given by

.. math::

   0 \le \rho \le \rho_\mathrm{max} \\
   \epsilon_c(\rho) \le \epsilon \le \epsilon_\mathrm{max} \\
   0 \le Y_e \le 1

where :math:`\rho_\mathrm{max}` and :math:`\epsilon_\mathrm{max}`
are free parameters specified by the user.

Since the thermal contribution in :math:numref:`def_hybrideos` is an 
ad-hoc prescription not derived from any nuclear physics model,
temperature and entropy are not uniquely defined (although one
can construct both by prescribing a functional dependence of
entropy and energy at fixed density). In this version of the library,
temperature and entropy are not implemented.



