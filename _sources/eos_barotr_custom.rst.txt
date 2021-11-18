Custom EOS
==========

It is easy to implement a custom EOS without modifying the
RePrimAnd library. One simply needs to define a class which
provides the EOS functions. This class needs to be derived
from a generic implementation interface
:cpp:class:`~EOS_Toolkit::implementations::eos_barotr_impl` 
and override the virtual functions defined by the interface.

Mostly, the implementation needs to be provide functions 
in terms of the pseudo-enthalpy :math:`g`, as well as a 
conversion between :math:`g` and :math:`\rho`:

.. math::

   &g(\rho) \\
   &\rho(g)  \\
   &P(g) \\
   &\epsilon(g) \\
   &h(g) \\
   &c_s(\rho, \theta, Y_e) \\
   &T(g) &\mathrm{(optional)} \\
   &Y_e(g) &\mathrm{(optional)} \\
   
In addition, it has to provide the validity range for
:math:`\rho` and :math:`g`, 
and the global minimum :math:`h_0` of the enthalpy.
If temperature and/or electron fraction are not provided, the 
corresponding methods should throw an exception. **Under no 
circumstances** should incorrect "dummy" values or NANs be returned. 


Here is an example declaration of a custom EOS implementation

.. code-block:: cpp
   
   #ifndef EOS_BARFOO_IMPL_H
   #define EOS_BARFOO_IMPL_H
   
   #include "eos_barotropic_impl.h"
   
   namespace EOS_Toolkit {
   namespace implementations {
   
   class eos_barfoo : public eos_barotr_impl {
     range rgrho;
     range rggm1;
     const real_t min_h;
     
   
     public:
   
     eos_barfoo(/* eos parameters */);
   
     const range& range_rho() const final {return rgrho;}
     
     const range& range_gm1() const final {return rggm1;}
   
     real_t minimal_h() const final {return min_h;}
   
     bool is_isentropic() const final;
     
     bool is_zero_temp() const final;
     
     bool has_temp() const final;
     
     bool has_efrac() const final;
   
   
     real_t gm1_from_rho(real_t rho) const final;
   
     real_t rho(real_t gm1) const final;
   
     real_t eps(real_t gm1) const  final; 
   
     real_t press(real_t gm1) const final;
   
     real_t hm1(real_t gm1) const final;
   
     real_t csnd(real_t gm1) const final;
   
     real_t temp(real_t gm1) const final;
   
     real_t ye(real_t gm1) const final;
   
   };
   
   }
   }
   
   #endif
   
 

No function will ever be called outside the validity range, all checks
are taken care of by the user-facing interface 
:cpp:class:`~EOS_Toolkit::eos_barotr`.
Conversely, the implementation should always return a correct result
for valid input, including parameters on the boundary of the valid 
region. 


Finally, one has to provide a function to wrap the
implementation into :cpp:class:`~EOS_Toolkit::eos_barotr` EOS object,
like this:


.. code-block:: cpp

   eos_barotr EOS_Toolkit::make_eos_barfoo(/* <barfoo parameters> */)
   {
     return eos_barotr{ 
             make_shared<eos_barfoo>(/* <barfoo parameters> */) };
   }

The custom EOS is then ready to use:

.. code-block:: cpp

   auto barfoo = make_eos_barfoo(/* <barfoo parameters> */);   
   press = barfoo.at_rho(rho).press();
   
   
For a real example, we suggest to look at the implementation of 
the polytropic EOS (under `library/EOS_Barotropic`)

Extending EOS file format
^^^^^^^^^^^^^^^^^^^^^^^^^

The library provides a mechanism to register a file reader for 
custom EOS, without changing the library itself. 
The universal EOS file format is open in the sense that
all EOS-type specific information is contained in a HDF5 subgroup
and the file has a string attribute `eos_type` for the type of the EOS.

When creating a file for a custom EOS `foobar`, the eos type should be 
named `barotr_custom_foobar` and the group holding EOS data should be 
named `eos_barotr_custom_foobar`. 

To register a reader, one needs to create a translation unit similar
to the one below. Files with custom EOS can then be loaded via the same
interface as for the types provided by the library.

.. code:: cpp

   #include "hdf5imple.h"
   #include "eos_barotr_file_impl.h"
   #include "eos_barotr_poly.h"
   
   namespace EOS_Toolkit {
   namespace implementations {
   
   
   struct reader_eos_barfoo : reader_eos_barotr 
   {
     eos_barotr load(const h5grp& g, const units& u) const final;
   };
   
   const bool register_reader_eos_barfoo { 
     registry_reader_eos_barotr::add("barotr_custom_barfoo", 
                                     new reader_eos_barfoo())
   };
   
   eos_barotr reader_eos_barfoo::load(const h5grp& g, 
                                           const units& u) const
   {
     
     // code to read EOS data from HDF5 group g goes here
     
     return make_eos_barfoo(/* barfoo EOS parameters loaded above */);
   }    
     
   } 
   }
   

The header `hdf5imple.h` provides a minimalistic C++ wrapper of the 
HDF5 interface, but one can also use hdf5 directly. To get the hdf5
handle of the group g, use `g.use()`. The file readers for existing
EOS are implemented in the same way as above and may serve as examples.



Reference
^^^^^^^^^

.. doxygenclass:: EOS_Toolkit::implementations::eos_barotr_impl
   :project: RePrimAnd
   :members:
