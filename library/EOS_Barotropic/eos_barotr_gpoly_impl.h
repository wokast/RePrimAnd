#ifndef EOS_BAROTR_GPOLY_IMPL_H
#define EOS_BAROTR_GPOLY_IMPL_H

#include "eos_barotropic_impl.h"

namespace EOS_Toolkit {
namespace implementations {

///Generalized Polytropic EOS
/**
The generalized polytropic EOS is the same as a polytropic one, except
that the specific internal energy is not assumed to approach zero are
zero density, but a value \f$ \epsilon_0 \f$. This allows the use 
of arbitrary conventions for the baryon mass.

See eos_cold for notation used and eos_cold_api for a description of 
the member functions.
*/
class eos_barotr_gpoly : public eos_barotr_impl {
  range rgrho;
  range rggm1;
  
  real_t n;       ///< Polytropic index \f$ n \f$
  real_t rmd_p;   ///< Polytropic density scale \f$ \rho_p \f$
  real_t np1;     ///< \f$ n+1 \f$
  real_t gamma;   ///< Polytropic exponent \f$ \Gamma \f$
  real_t invn;    ///< \f$ \frac{1}{n} \f$
  real_t sed0;    ///< \f$ \epsilon_0 \f$
  real_t h0;      ///< \f$ h_0 = 1 + \epsilon_0 \f$


  public:

///Constructor
  eos_barotr_gpoly(
    real_t n_,                          ///<Adiabatic index \f$ n \f$
    real_t rmd_p_,                      ///<Density scale \f$ \rho_p \f$
    real_t sed0_,                       ///< \f$ \epsilon_0 \f$    
    real_t rho_max_                     ///<Max valid density 
  );
  
  eos_barotr_gpoly(const eos_barotr_gpoly&)            = default;
  eos_barotr_gpoly(eos_barotr_gpoly&&)                 = default;
  ~eos_barotr_gpoly() final                            = default;
  eos_barotr_gpoly& operator=(eos_barotr_gpoly&&)      = delete;
  eos_barotr_gpoly& operator=(const eos_barotr_gpoly&) = delete;


  ///Returns range of validity for density
  const range& range_rho() const final {return rgrho;}
  
  ///Returns range of validity for \f$ g-1 \f$ 
  const range& range_gm1() const final {return rggm1;}

  ///Returns range of validity for \f$ g-1 \f$ 
  real_t minimal_h() const final {return h0;}

  ///Whether EOS is isentropic
  bool is_isentropic() const final  {return true;}
  
  ///Whether EOS is for zero temperature
  bool is_zero_temp() const final  {return true;}
  
  ///Whether EOS can compute temperature
  bool has_temp() const final  {return true;}
  
  ///Whether EOS can compute electron fraction
  bool has_efrac() const final  {return false;}  


  ///Compute \f$ g-1 \f$ 
  /**Assumes input is in the valid range, no checks are performed.*/
  real_t gm1_from_rho(
    real_t rho      ///<Rest mass density  \f$ \rho \f$
  ) const final;

  ///Compute Rest mass density \f$ \rho \f$
  /**Assumes input is in the valid range, no checks are performed.*/
  real_t rho(
    real_t gm1      ///< \f$ g-1 \f$
  ) const final;


  ///Compute Specific internal energy \f$\epsilon \f$
  /**Assumes input is in the valid range, no checks are performed.*/
  real_t eps(
    real_t gm1      ///< \f$ g-1 \f$
  ) const final; 


  ///Compute Pressure \f$ P \f$ 
  /**Assumes input is in the valid range, no checks are performed.*/
  real_t press(
    real_t gm1      ///< \f$ g-1 \f$
  ) const final;


  ///Compute specific enthalpy (excluding restmass) \f$ h-1  \f$
  /**Assumes input is in the valid range, no checks are performed.*/
  real_t hm1(
    real_t gm1      ///< \f$ g-1 \f$
  ) const final;

  ///Compute adiabatic soundspeed \f$ c_s \f$
  /**Assumes input is in the valid range, no checks are performed.*/
  real_t csnd(
    real_t gm1      ///< \f$ g-1 \f$
  ) const final;

  ///Returns temperature \f$ T = 0\f$ for this EOS 
  /**Assumes input is in the valid range, no checks are performed.*/
  real_t temp(
    real_t gm1      ///< \f$ g-1 \f$
  ) const final {return 0.0;}

  ///Computing electron fraction not implemented for this EOS
  real_t ye(
    real_t gm1      ///< \f$ g-1 \f$
  ) const final;

  static real_t rmd_p_from_p_rho_n(real_t p, real_t rho, real_t n);
  static real_t eps0_from_p_rho_eps_n(real_t p, real_t rho, 
                                      real_t eps, real_t n);

  static eos_barotr_gpoly 
  from_boundary(real_t rho0, real_t eps0, real_t p0, 
                real_t n_poly, real_t rho_max);
 
};

}// namespace implementations 
}// namespace EOS_Toolkit 


#endif

