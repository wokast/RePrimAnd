
Reference
---------

Representation of Sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenclass:: EOS_Toolkit::star_seq
   :project: RePrimAnd
   :members:

|


.. doxygenclass:: EOS_Toolkit::star_branch
   :project: RePrimAnd
   :members:

Creating star sequences
^^^^^^^^^^^^^^^^^^^^^^^


Computing TOV Sequences
~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: EOS_Toolkit::make_tov_seq(eos_barotr eos, interval<real_t> rg_gm1, const star_accuracy_spec acc, unsigned int num_samp)
   :project: RePrimAnd

|

.. doxygenfunction:: EOS_Toolkit::make_tov_branch_stable(eos_barotr eos, const star_accuracy_spec acc, real_t mg_cut_low_rel=0.2, real_t mg_cut_low_abs=0.0, real_t gm1_initial=1.2, real_t gm1_step=0.004, real_t max_margin=1e-2)
   :project: RePrimAnd


 
Create from existing data
~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: EOS_Toolkit::make_star_seq(std::vector<real_t> mg, std::vector<real_t> mb, std::vector<real_t> rc, std::vector<real_t> mi, std::vector<real_t> lt, star_seq::range_t rg_gm1, units u)
   :project: RePrimAnd


Sequence Files 
^^^^^^^^^^^^^^

File operations
~~~~~~~~~~~~~~~

.. doxygenfunction:: EOS_Toolkit::save_star_seq(std::string fname, const star_seq& t)
   :project: RePrimAnd

|

.. doxygenfunction:: EOS_Toolkit::save_star_branch(std::string fname, const star_branch& t)
   :project: RePrimAnd

|

.. doxygenfunction:: EOS_Toolkit::load_star_seq(std::string fname, const units& u)
   :project: RePrimAnd

|

.. doxygenfunction:: EOS_Toolkit::load_star_branch(std::string fname, const units& u)
   :project: RePrimAnd


