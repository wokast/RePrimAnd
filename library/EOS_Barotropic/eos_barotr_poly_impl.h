#ifndef EOS_BAROTR_POLY_IMPL_H
#define EOS_BAROTR_POLY_IMPL_H

#include "eos_barotropic_impl.h"

namespace EOS_Toolkit {
namespace implementations {

///Polytropic EOS
/**
The polytropic EOS is written (using units with \f$ c=1 \f$) as
\f[ P = \rho_p \left(\frac{\rho}{\rho_p}\right)^\Gamma ,\qquad
\Gamma= 1+\frac{1}{n}  \f]
Here we use a <em>polytropic density scale</em> \f$ \rho_p \f$ to specify
the EOS instead of the usual form \f$ P = K \rho^\Gamma \f$
because it has simpler units than 
\f$ K = \rho_p^{-1/n} \f$.
*/
class eos_barotr_poly : public eos_barotr_impl {
  range rgrho;
  range rggm1;
  const real_t min_h{1.0};
  
  real_t n;       ///< Polytropic index \f$ n \f$
  real_t rmd_p;   ///< Polytropic density scale \f$ \rho_p \f$
  real_t np1;     ///< \f$ n+1 \f$
  real_t gamma;   ///< Polytropic exponent \f$ \Gamma \f$
  real_t invn;    ///< \f$ \frac{1}{n} \f$

  void init(
    real_t n_,                          ///<Adiabatic index \f$ n \f$
    real_t rmd_p_,                      ///<Density scale \f$ \rho_p \f$
    real_t rho_max_                     ///<Max valid density 
  );

  public:

  ///Constructor
  eos_barotr_poly(
    real_t n_,                          ///<Adiabatic index \f$ n \f$
    real_t rmd_p_,                      ///<Density scale \f$ \rho_p \f$
    real_t rho_max_                     ///<Max valid density 
  );

  ///Alternative constructor
  eos_barotr_poly(
    real_t rmd_m, ///<Rest mass density \f$ \rho \f$ at matching point
    real_t sed_m, ///<Specific internal energy \f$ \epsilon \f$ at matching point
    real_t p_m,   ///<Pressure \f$ P \f$ at matching point
    real_t rho_max_ ///<Max valid density 
  );


  ///Returns range of validity for density
  const range& range_rho() const final {return rgrho;}
  
  ///Returns range of validity for \f$ g-1 \f$ 
  const range& range_gm1() const final {return rggm1;}

  ///Returns range of validity for \f$ g-1 \f$ 
  real_t minimal_h() const final {return min_h;}

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
  ) const  final; 


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

};

}// namespace implementations
}//namespace EOS_Toolkit


#endif

