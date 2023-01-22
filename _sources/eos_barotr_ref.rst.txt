Reference
---------


Loading from File
^^^^^^^^^^^^^^^^^

.. doxygenfunction:: EOS_Toolkit::load_eos_barotr
   :project: RePrimAnd


Creating EOS from sampled data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To create a EOS based on spline interpolation from given sample points, 
there are three options.
Depending on whether the pseudo-enthalpy and/or specific energy are 
known accurately, one should use one of the following.

.. doxygenfunction:: EOS_Toolkit::make_eos_barotr_spline(const std::vector<real_t>& gm1, const std::vector<real_t>& rho, const std::vector<real_t>& eps, const std::vector<real_t>& press, const std::vector<real_t>& csnd, const std::vector<real_t>& temp, const std::vector<real_t>& efrac, bool isentropic, interval<real_t> rg_rho, real_t n_poly, units u, std::size_t pts_per_mag)
   :project: RePrimAnd

|

.. doxygenfunction:: EOS_Toolkit::make_eos_barotr_spline(const std::vector<real_t>& rho, const std::vector<real_t>& eps, const std::vector<real_t>& press, const std::vector<real_t>& csnd, const std::vector<real_t>& temp, const std::vector<real_t>& efrac, bool isentropic, interval<real_t> rg_rho, real_t n_poly, units u, std::size_t pts_per_mag)
   :project: RePrimAnd

|

.. doxygenfunction:: EOS_Toolkit::make_eos_barotr_spline(const std::vector<real_t>& rho, const std::vector<real_t>& press, const std::vector<real_t>& csnd, const std::vector<real_t>& temp, const std::vector<real_t>& efrac, interval<real_t> rg_rho, real_t n_poly, real_t eps0, units u, std::size_t pts_per_mag)
   :project: RePrimAnd

|

It is also possible to sample an existing EOS:

.. doxygenfunction:: EOS_Toolkit::make_eos_barotr_spline(const eos_barotr& eos, interval<real_t> rg_rho, real_t n_poly, std::size_t pts_per_mag=200)
   :project: RePrimAnd

|

Creating Polytropic EOS
^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: EOS_Toolkit::make_eos_barotr_poly
   :project: RePrimAnd

|

Creating Piecewise Polytropic EOS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: EOS_Toolkit::make_eos_barotr_pwpoly
   :project: RePrimAnd



Generic Interface
^^^^^^^^^^^^^^^^^
.. doxygenclass:: EOS_Toolkit::eos_barotr
   :project: RePrimAnd
   :members:|


Deprecated EOS 
^^^^^^^^^^^^^^

In older versions, there was no spline based EOS. Instead 
there was eos_barotr_table based on linear interpolation. This type 
can still be used for sake of reproducing old results, but for new
EOS its use is discouraged.

.. deprecated:: 1.5
   Use make_eos_barotr_spline instead make_eos_barotr_table

.. doxygenfunction:: EOS_Toolkit::make_eos_barotr_table
   :project: RePrimAnd
