
Reference
---------

Computing TOV solutions
^^^^^^^^^^^^^^^^^^^^^^^


Full solution
~~~~~~~~~~~~~

.. doxygenfunction:: EOS_Toolkit::make_tov_star(const eos_barotr eos, const real_t rho_center, const tov_acc_simple acc, const bool find_bulk, const bool find_tidal)
   :project: RePrimAnd

|

.. doxygenfunction:: EOS_Toolkit::make_tov_star(const eos_barotr eos, const real_t rho_center, const tov_acc_precise acc, const bool find_bulk, const bool find_tidal)
   :project: RePrimAnd

Only global properties 
~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: EOS_Toolkit::get_tov_star_properties(const eos_barotr eos, const real_t rho_center, const tov_acc_simple acc, const bool find_bulk, const bool find_tidal)
   :project: RePrimAnd

|

.. doxygenfunction:: EOS_Toolkit::get_tov_star_properties(const eos_barotr eos, const real_t rho_center, const tov_acc_precise acc, const bool find_bulk, const bool find_tidal)
   :project: RePrimAnd


NS representations
^^^^^^^^^^^^^^^^^^

.. doxygenclass:: EOS_Toolkit::spherical_star_properties
   :project: RePrimAnd
   :members:

|

.. doxygenclass:: EOS_Toolkit::spherical_star
   :project: RePrimAnd
   :members:


Accuracy specification
^^^^^^^^^^^^^^^^^^^^^^

.. doxygenstruct:: EOS_Toolkit::tov_acc_simple
   :project: RePrimAnd
   :members:

|

.. doxygenstruct:: EOS_Toolkit::tov_acc_precise
   :project: RePrimAnd
   :members:



TOV sequences
^^^^^^^^^^^^^

.. doxygenfunction:: EOS_Toolkit::find_rhoc_tov_max_mass(eos_barotr eos, const real_t rhobr0, const real_t rhobr1, const int bits, const real_t acc, unsigned int max_steps)
   :project: RePrimAnd

|

.. doxygenfunction:: EOS_Toolkit::find_rhoc_tov_of_mass(eos_barotr eos, real_t mg, const real_t rhobr0, const real_t rhobr1, real_t acc, unsigned int max_steps)
   :project: RePrimAnd


