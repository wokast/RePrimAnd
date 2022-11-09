Reference
---------


Generic Interface
^^^^^^^^^^^^^^^^^
.. doxygenclass:: EOS_Toolkit::eos_barotr
   :project: RePrimAnd
   :members:

Loading from File
^^^^^^^^^^^^^^^^^

.. doxygenfunction:: EOS_Toolkit::load_eos_barotr
   :project: RePrimAnd


Creating Specific EOS
^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: EOS_Toolkit::make_eos_barotr_spline(const std::vector<real_t>& gm1, const std::vector<real_t>& rho, const std::vector<real_t>& eps, const std::vector<real_t>& press, const std::vector<real_t>& csnd, const std::vector<real_t>& temp, const std::vector<real_t>& efrac, bool isentropic_, interval<real_t> rg_rho, real_t n_poly_, units units_, std::size_t pts_per_mag)
   :project: RePrimAnd

|

.. doxygenfunction:: EOS_Toolkit::make_eos_barotr_spline(const eos_barotr& eos, interval<real_t> rg_rho, real_t n_poly, std::size_t pts_per_mag=200)
   :project: RePrimAnd

|

.. doxygenfunction:: EOS_Toolkit::make_eos_barotr_poly
   :project: RePrimAnd

|

.. doxygenfunction:: EOS_Toolkit::make_eos_barotr_pwpoly
   :project: RePrimAnd

|


.. deprecated:: 1.5
   Use make_eos_barotr_spline instead make_eos_barotr_table

.. doxygenfunction:: EOS_Toolkit::make_eos_barotr_table
   :project: RePrimAnd
