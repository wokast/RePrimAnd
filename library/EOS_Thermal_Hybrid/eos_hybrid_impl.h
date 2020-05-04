#ifndef EOSHYBRID_IMPL_H
#define EOSHYBRID_IMPL_H
#include "eos_barotropic.h"
#include "eos_thermal_impl.h"

namespace EOS_Toolkit {
namespace implementations {

class eos_hybrid : public eos_thermal_impl {
  eos_barotr eos_c;
  real_t gamma_th, gm1_th, eps_max;
  range rgrho;      ///< Valid range for density \f$ \rho \f$
  range rgye;       ///< Valid range for electron fraction \f$ Y_e \f$
  real_t min_h;    ///< Lower bound for enthalpy \f$ h \ge h_0 > 0 \f$


  real_t eps_cold(real_t rho) const;
  real_t p_cold(real_t rho) const;
  real_t hm1_cold(real_t rho) const;
  real_t cs2_cold(real_t rho) const;

  public:

  eos_hybrid(eos_barotr eos_c_,  real_t gamma_th_,  
             real_t eps_max_, real_t rho_max_); 
  ~eos_hybrid() final = default;


  ///Identity: thermal variable is eps for this EOS
  real_t therm_from_rho_eps_ye(
    real_t rho,     ///<Rest mass density  \f$ \rho \f$
    real_t eps,     ///Specific internal energy \f$ \epsilon \f$
    real_t ye       ///<Electron fraction \f$ Y_e \f$ 
  ) const final
  {
    return eps;
  }

  ///Compute specific energy from temperature.
  real_t therm_from_rho_temp_ye(
    real_t rho,     ///<Rest mass density  \f$ \rho \f$
    real_t temp,    ///Temperature
    real_t ye       ///<Electron fraction \f$ Y_e \f$ 
  ) const final;

  ///Inverse Identity: eps is thermal variable for this EOS
  real_t eps(
    real_t rho,     ///<Rest mass density  \f$ \rho \f$
    real_t eps,     ///<Specific internal energy \f$ \epsilon \f$
    real_t ye       ///<Electron fraction \f$ Y_e \f$ 
  ) const final 
  {
    return eps;
  }
  
  ///Compute temperature
  real_t temp(
    real_t rho,     ///<Rest mass density  \f$ \rho \f$
    real_t eps,     ///<Specific internal energy \f$ \epsilon \f$
    real_t ye       ///<Electron fraction \f$ Y_e \f$ 
  ) const final;

  ///Compute pressure
  real_t press(
    real_t rho,     ///<Rest mass density  \f$ \rho \f$
    real_t eps,     ///<Specific internal energy \f$ \epsilon \f$
    real_t ye       ///<Electron fraction \f$ Y_e \f$ 
  ) const final;

  ///Compute soundspeed
  real_t csnd(
    real_t rho,     ///<Rest mass density  \f$ \rho \f$
    real_t eps,     ///<Specific internal energy \f$ \epsilon \f$
    real_t ye       ///<Electron fraction \f$ Y_e \f$ 
  ) const final;

  ///Compute entropy. Not implemented.
  real_t sentr(
    real_t rho,     ///<Rest mass density  \f$ \rho \f$
    real_t eps,     ///<Specific internal energy \f$ \epsilon \f$
    real_t ye       ///<Electron fraction \f$ Y_e \f$ 
  ) const final;

  real_t dpress_drho(
    real_t rho,     ///<Rest mass density  \f$ \rho \f$
    real_t eps,     ///<Specific internal energy \f$ \epsilon \f$
    real_t ye       ///<Electron fraction \f$ Y_e \f$ 
  ) const final;

  real_t dpress_deps(
    real_t rho,     ///<Rest mass density  \f$ \rho \f$
    real_t eps,     ///<Specific internal energy \f$ \epsilon \f$
    real_t ye       ///<Electron fraction \f$ Y_e \f$ 
  ) const final;


  /// Valid range for density
  const range& range_rho() const final {return rgrho;}

  /// Valid range for electron fraction
  const range& range_ye() const final {return rgye;}

  range range_eps(
    real_t rho,     ///<Rest mass density  \f$ \rho \f$
    real_t ye       ///<Electron fraction \f$ Y_e \f$ 
  ) const final;

  range range_temp(
    real_t rho,     ///<Rest mass density  \f$ \rho \f$
    real_t ye       ///<Electron fraction \f$ Y_e \f$ 
  ) const final;

  real_t minimal_h() const final {return min_h;}

};

} // namespace implementations

} // namespace EOS_Toolkit

#endif


