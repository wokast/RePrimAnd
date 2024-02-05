
Reference
---------

Computing TOV solutions
^^^^^^^^^^^^^^^^^^^^^^^


Full solution
~~~~~~~~~~~~~

.. doxygenfunction:: EOS_Toolkit::get_tov_star(const eos_barotr eos, const real_t rho_center, const star_accuracy_spec acc)
   :project: RePrimAnd

|

Only global properties 
~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: EOS_Toolkit::get_tov_properties(const eos_barotr eos, const real_t rho_center, const star_accuracy_spec acc)
   :project: RePrimAnd

|


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

.. doxygenstruct:: EOS_Toolkit::star_accuracy_spec
   :project: RePrimAnd
   :members:

|

.. doxygenfunction:: EOS_Toolkit::star_acc_simple(bool need_deform, bool need_bulk,  real_t acc_tov, real_t acc_deform, std::size_t minsteps)
   :project: RePrimAnd

.. doxygenfunction:: EOS_Toolkit::star_acc_detailed(bool need_deform, bool need_bulk,  real_t acc_mass, real_t acc_radius, real_t acc_minertia, real_t acc_deform, std::size_t minsteps) 
   :project: RePrimAnd




Other
^^^^^

The following low level functions should rarely be needed, usually
it is more convenient to use the star sequence functionality.

.. doxygenfunction:: EOS_Toolkit::find_rhoc_tov_max_mass(eos_barotr eos, const real_t rhobr0, const real_t rhobr1, const int bits, const real_t acc, unsigned int max_steps)
   :project: RePrimAnd

|

.. doxygenfunction:: EOS_Toolkit::find_rhoc_tov_of_mass(eos_barotr eos, real_t mg, const real_t rhobr0, const real_t rhobr1, real_t acc, unsigned int max_steps)
   :project: RePrimAnd


