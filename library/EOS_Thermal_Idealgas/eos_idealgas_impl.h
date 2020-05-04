#ifndef EOS_IDEALGAS_IMPL_H
#define EOS_IDEALGAS_IMPL_H

#include "eos_thermal_impl.h"

namespace EOS_Toolkit {

namespace implementations {

///This implements an ideal gas EOS.
/**
\note
The EOS framework includes electron fraction, but here it is ignored 
and serves as a dummy parameter. 
The temperature is the one from the classic ideal gas, not some 
re-interpretation where it is zero along some polytrope. The entropy is 
currently not implemented.
**/
class eos_idealgas : public eos_thermal_impl {
  real_t gamma, gm1;
  range rgrho;      ///< Valid range for density \f$ \rho \f$
  range rgye;       ///< Valid range for electron fraction \f$ Y_e \f$
  real_t min_h;    ///< Lower bound for enthalpy \f$ h \ge h_0 > 0 \f$

  range rgeps;
  
  public:

  eos_idealgas(real_t n_, real_t max_eps_, real_t max_rho_);
                
  ~eos_idealgas() final = default;

  ///Identity: thermal variable is eps for this EOS
  real_t therm_from_rho_eps_ye(
    real_t rho,     ///<Rest mass density  \f$ \rho \f$
    real_t eps,     ///<Specific internal energy \f$ \epsilon \f$
    real_t ye       ///<Electron fraction \f$ Y_e \f$ 
  ) const final 
  {
    return eps;
  }

  ///Compute specific energy from temperature.
  [[ noreturn ]] 
  real_t therm_from_rho_temp_ye(
    real_t rho,     ///<Rest mass density  \f$ \rho \f$
    real_t temp,    ///<Temperature
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
  
  ///Compute temperature. Not implemented.
  [[ noreturn ]] 
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
  [[ noreturn ]] 
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

  [[ noreturn ]] 
  range range_temp(
    real_t rho,     ///<Rest mass density  \f$ \rho \f$
    real_t ye       ///<Electron fraction \f$ Y_e \f$ 
  ) const final;

  real_t minimal_h() const final {return min_h;}

};


} // namespace implementations
} // namespace EOS_Toolkit

#endif


