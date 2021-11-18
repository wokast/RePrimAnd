Notation
========

Units
^^^^^
Unless specifically noted, we use geometric units :math:`G=c=1`,
both for documentation and the code internals. Which specific 
geometric units to employ is up to the user, e.g. setting 
:math:`M_\odot=1`. 
All EOS and primitive recovery code is completely
agnostic to this choice. 
Only when loading an EOS from file, the unit system needs to be 
specified. Our EOS file format is using SI units in order to 
unambiguously define an EOS.

We further define the fiducial baryon 
mass constant (used to convert baryon number density to mass density) 
as :math:`m_B = 1.66 \times 10^{-24}\, \mathrm{g}`.
This is mostly a convention to make equations of state unambiguous, 
which is neccessary because our EOS interface is based on baryon mass 
density, not at all using baryon number density.

Temperatures are always measured in :math:`[\mathrm{MeV}]` and specific 
entropy in :math:`[k_B/\mathrm{Baryon}]`.
Computing temperature/entropy is the only place where the value of 
:math:`m_B` might be used in the code.


Primitive variables
^^^^^^^^^^^^^^^^^^^

.. list-table:: 
   :widths: 30 20 50
   :header-rows: 1

   * - Math
     - Code
     - Explanation
   * - :math:`n`
     -
     - Baryon number density
   * - :math:`m_B`
     - ``mbar``
     - Fiducial baryon mass constant
   * - :math:`\rho = m_B n`
     - ``rho``
     - Mass density
   * - :math:`\rho_E`
     - 
     - Fluid contribution to energy density
   * - :math:`\epsilon=\frac{\rho_E}{\rho} - 1`
     - ``eps`` 
     - Fluid specific internal energy
   * - :math:`P`
     - ``press``
     - Fluid pressure
   * - :math:`h=1+\epsilon + \frac{P}{\rho}`
     - ``h``
     - Relativistic specific enthalpy
   * - :math:`T`
     - ``temp``
     - Temperature
   * - :math:`s`
     - ``sentr``
     - Entropy per baryon
   * - :math:`c_s = \sqrt{\left.\frac{\partial P}{\partial \rho_E} \right|_s}`
     - ``csnd``
     - Speed of sound (adiabatic)
   * - :math:`v^i`
     - ``vel``
     - Fluid 3-velocity for Eularian frame (observer normal to foliation)
   * - :math:`v^2 = v^i v_i`
     - ``vsqr``
     -
   * - :math:`W`
     - ``w_lor``
     - Fluid Lorentz factor in Eularian frame
   * - :math:`z = Wv`
     - ``z``
     - More useful quantity for expressing velocity
   * - :math:`B^i`
     - ``B``
     - Magnetic field in Eularian frame
   * - :math:`E^i`
     - ``E``
     - Electric field in Eularian frame

Magnetic and electric field are defined as 

.. math::
   E^\mu &= n_\nu F^{\mu\nu}, \\
   B^\mu &= n_\nu {}^*F^{\mu\nu}

where :math:`F^{\mu\nu}` is the Maxwell tensor and the star denotes the dual. Beware of 
competing conventions in the literature (rationalized or not). The energy density 
is given by

.. math::
   \frac{1}{2} \left( E^2 + B^2 \right)

Evolved Variables
^^^^^^^^^^^^^^^^^
The definition of the evolved (a.k.a. conserved) variables follows 
typical evolution code conventions. All evolved variables are tensor 
densities, i.e. incorporate a factor :math:`\sqrt{g}`, the square root 
of the 3-metric determinant. We also use that convention for the 
evolved magnetic field. Further, we assume that the electron fraction
is evolved as a densitized tracer variable. Conversion from other 
conventions should be trivial.

.. list-table::
   :widths: 50 20 30
   :header-rows: 1

   * - Math
     - Code
     - Explanation
   * - :math:`D = \sqrt{g} \rho W`
     - ``dens``
     - Evolved mass density
   * - :math:`\tau = \sqrt{g}\left( \rho W \left(hW-1\right) - P + \frac{1}{2} \left(E^2 + B^2 \right) \right)`
     - ``tau``
     - Evolved energy density
   * - :math:`S_i = \sqrt{g}\left( \rho h W^2 v_i + \epsilon_{ijk} E^j B^k \right)`
     - ``scon``
     - Evolved momentum density 
   * - :math:`Y_e^T = D Y_e`
     - ``tracer_ye``
     - Evolved electron fraction tracer
   * - :math:`B_c^i = \sqrt{g} B^i`
     - ``bcons``
     - Evolved magnetic field 
     
