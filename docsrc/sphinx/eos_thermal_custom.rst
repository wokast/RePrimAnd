Custom EOS
==========

It is easy to implement a custom EOS without modifying the
RePrimAnd library. One simply needs to define a class which
provides the EOS functions. This class needs to be derived
from a generic implementation interface
:cpp:class:`~EOS_Toolkit::implementations::eos_thermal_impl` 
and override the virtual functions defined by the interface.

One subtlety regards the representation of the thermal degree
of freedom. For some EOS, it might be more natural to use
temperature, for others specific energy or something else entirely. 
In the following, we denote the "natural" choice as :math:`\theta`.
We leave the choice of :math:`\theta` open, as an implementation 
detail. 

The implementation needs to provide the following functions

.. math::

   &\theta(\rho, \epsilon, Y_e) \\
   &\theta(\rho, T, Y_e) &\mathrm{(optional)} \\
   &P(\rho, \theta, Y_e) \\
   &\epsilon(\rho, \theta, Y_e) \\
   &T(\rho,\theta, Y_e) &\mathrm{(optional)} \\
   &s(\rho,\theta, Y_e) &\mathrm{(optional)} \\
   &c_s(\rho, \theta, Y_e) \\
   &\frac{\partial P}{\partial \rho}(\rho, \theta, Y_e) \\
   &\frac{\partial P}{\partial \epsilon}(\rho, \theta, Y_e)

In addition, it has to provide the validity range as
:math:`\rho_\mathrm{min}, \rho_\mathrm{max}`, 
:math:`Y_{e,\mathrm{min}}, Y_{e,\mathrm{max}}`, 
:math:`\epsilon_\mathrm{min}(\rho, Y_e), \epsilon_\mathrm{max}(\rho, Y_e)`,
and the global minimum :math:`h_0` of the enthalpy.

Here is an example declaration of a custom EOS implementation

.. code-block:: cpp

   #ifndef EOS_FOOBAR_IMPL_H
   #define EOS_FOOBAR_IMPL_H
 
   #include "eos_thermal_impl.h"
 
   namespace EOS_Toolkit {
 
   namespace implementations {
     
   class eos_foobar : public eos_thermal_impl {
     range rgrho;      
     range rgye;       
     real_t min_h; 
 

     public:
 
     eos_foobar(/* foobar parameters */);
                   
     virtual ~eos_foobar();
 
     real_t therm_from_rho_eps_ye(real_t rho,real_t eps, 
                                  real_t ye) const override; 
 
     real_t therm_from_rho_temp_ye(real_t rho, real_t temp, 
                                   real_t ye) const override;
 
     real_t eps(real_t rho, real_t therm, 
                real_t ye) const override;
     
     real_t temp(real_t rho, real_t therm, 
                 real_t ye) const override;
 
     real_t press(real_t rho, real_t therm, 
                  real_t ye) const override;
 
     real_t csnd(real_t rho, real_t therm, 
                 real_t ye) const override;
 
     real_t sentr(real_t rho, real_t therm, 
                  real_t ye) const override;
 
     real_t dpress_drho(real_t rho, real_t therm, 
                        real_t ye) const override;
 
     real_t dpress_deps(real_t rho, real_t therm, 
                        real_t ye) const override;
     
     const range& range_rho() const override {return rgrho;}
     const range& range_ye() const override {return rgye;}
 
     range range_eps(real_t rho, real_t ye) const override;
     range range_temp(real_t rho, real_t ye) const override;
 
     real_t minimal_h() const override {return min_h;}
 
   };
 
 
   } 
   } 
   
   #endif
 

No function will ever be called outside the validity range, all checks
are taken care of by the user-facing interface 
:cpp:class:`~EOS_Toolkit::eos_thermal`.
Conversely, the implementation should always return a correct result
for valid input, including parameters on the boundary of the valid 
region. 

If temperature and/or entropy are not provided, the corresponding 
methods should throw an exception. **Under no circumstances** should
incorrect "dummy" values or NANs be returned. 

Finally, one has to provide a function to wrap the
implementation into :cpp:class:`~EOS_Toolkit::eos_thermal` EOS object,
like this:


.. code-block:: cpp

   eos_thermal EOS_Toolkit::make_eos_foobar(/* <foobar parameters> */)
   {
     return eos_thermal{ 
               make_shared<eos_foobar>(/* <foobar parameters> */) };
   }

The custom EOS is then ready to use:

.. code-block:: cpp

   auto foobar = make_eos_foobar(/* <foobar parameters> */);   
   press = foobar.at_rho_eps_ye(rho,eps,ye).press();
   
   
For a real example, we suggest to look at the implementation of 
the ideal gas EOS (under `library/EOS_Thermal_Idealgas`)

Extending EOS file format
^^^^^^^^^^^^^^^^^^^^^^^^^

The library provides a mechanism to register a file reader for 
custom EOS, without changing the library itself. 
The universal EOS file format is open in the sense that
all EOS-type specific information is contained in a HDF5 subgroup
and the file has a string attribute `eos_type` for the type of the EOS.

When creating a file for a custom EOS `foobar`, the eos type should be 
named `thermal_custom_foobar` and the group holding EOS data should be 
named `eos_thermal_custom_foobar`. 

To register a reader, one needs to create a translation unit similar
to the one below. Files with custom EOS can then be loaded via the same
interface as for the types provided by the library.

.. code:: cpp

   #include "hdf5imple.h"
   #include "eos_thermal_file_impl.h"
   #include "eos_foobar.h"
   
   namespace EOS_Toolkit {
   namespace implementations {
   
   
   struct reader_eos_foobar : reader_eos_thermal 
   {
     eos_thermal load(const h5grp& g, const units& u) const final;
   };
   
   const bool register_reader_eos_foobar { 
     registry_reader_eos_thermal::add("thermal_custom_foobar", 
                                      new reader_eos_foobar())
   };
   
   eos_thermal reader_eos_foobar::load(const h5grp& g, 
                                              const units& u) const
   {
     //code to read EOS from HDF5 group g goes here
                       
     return make_eos_foobar(/* foobar parameters loaded above */);
   }
     
     
     
   } //namespace implementations
   } //namespace EOS_Toolkit
   

The header `hdf5imple.h` provides a minimalistic C++ wrapper of the 
HDF5 interface, but one can also use hdf5 directly. To get the hdf5
handle of the group g, use `g.use()`. The file readers for existing
EOS are implemented in the same way as above and may serve as examples.


Reference
^^^^^^^^^

.. doxygenclass:: EOS_Toolkit::implementations::eos_thermal_impl
   :project: RePrimAnd
   :members:
