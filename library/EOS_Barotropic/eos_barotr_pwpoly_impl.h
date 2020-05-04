#ifndef EOS_BAROTR_PWPOLY_IMPL_H
#define EOS_BAROTR_PWPOLY_IMPL_H

#include "eos_barotropic_impl.h"
#include <vector>

namespace EOS_Toolkit {
namespace implementations {

///Class representing a segment of the piecewise polytropic EOS
class eos_poly_piece {
  public:
  real_t rmd0;
  real_t dsed;
  real_t gamma;
  real_t rmd_p;
  real_t n;
  real_t np1;
  real_t invn;
  real_t gm10;
  real_t p0;

  eos_poly_piece() = default;  
  eos_poly_piece(real_t rmd0_, real_t sed0_, 
                 real_t gamma_, real_t rmd_p_);



  real_t gm1_from_rho(real_t rho) const;
  real_t eps_from_gm1(real_t gm1) const;
  real_t eps_from_rho(real_t rho) const;
  real_t press_from_gm1(real_t gm1) const;
  real_t rho_from_gm1(real_t gm1) const;
  real_t hm1_from_gm1(real_t gm1) const;
  real_t csnd_from_gm1(real_t gm1) const;
  
  real_t rho_max_save(real_t rho_max) const;
  bool rho_save_up_to(real_t rho) const;
};



///Piecewise Polytropic EOS
class eos_barotr_pwpoly : public eos_barotr_impl {
  range rgrho;
  range rggm1;
  const real_t min_h{1.0};

  ///The polytropic segments
  std::vector<eos_poly_piece> segments;
  
  ///Find the segment responsible for a given mass density
  const eos_poly_piece& segment_for_rho(real_t rho) const;
  ///Find the segment responsible for a given \f$ g - 1\f$
  const eos_poly_piece& segment_for_gm1(real_t gm1) const;
  
  public:

  ///Constructor
  eos_barotr_pwpoly(
    real_t rmdp0,                           ///<First segment polytropic density scale
    const std::vector<real_t>& segm_bound,  ///<Densities of segment boundaries
    const std::vector<real_t>& segm_gamma,  ///<Segment gammas
    real_t rho_max_                         ///<EOS max valid density
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

};

}//namespace implementations
}//namespace EOS_Toolkit


#endif

