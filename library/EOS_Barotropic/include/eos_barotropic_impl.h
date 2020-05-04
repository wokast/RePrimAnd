#ifndef EOS_BAROTROPIC_IMPL_H
#define EOS_BAROTROPIC_IMPL_H

#include "config.h"
#include "intervals.h"

namespace EOS_Toolkit {
namespace implementations {


///Abstract interface for barotropic equation of state implementations.
class eos_barotr_impl {
  public:
  using range  = interval<real_t>;

  eos_barotr_impl()          = default;
  eos_barotr_impl(const eos_barotr_impl&) = default;
  eos_barotr_impl(eos_barotr_impl&&) = default;
  eos_barotr_impl& operator=(const eos_barotr_impl&) = delete;
  eos_barotr_impl& operator=(eos_barotr_impl&&) = delete;
  virtual ~eos_barotr_impl() = default;

  /**
  @return Whether EOS is isentropic
  **/
  virtual bool is_isentropic() const =0;
  
  /**
  @return Whether EOS is for zero temperature
  **/
  virtual bool is_zero_temp() const =0;
  
  /**
  @return Whether EOS can compute temperature
  **/
  virtual bool has_temp() const =0;
  
  /**
  @return Whether EOS can compute electron fraction
  **/
  virtual bool has_efrac() const =0;
  
  /**
  @return Range of validity for mass density \f$ \rho \f$
  **/
  virtual const range& range_rho() const =0;
  
  /**
  @return Range of validity for \f$ g-1 \f$ 
  **/
  virtual const range& range_gm1() const =0;

  /**
  @return Global minimum \f$ h_0 \f$ of enthalpy.
  
  It must be guaranteed that the enthaply computed by the EOS for 
  any valid state is always above this value. 
  
  \post The minimum must be strictly positive.
  **/
  virtual real_t minimal_h() const =0;

  /**
  @param rho Rest mass density  \f$ \rho \f$
  @return Pseudo enthalpy \f$ g-1 \f$
  **/ 
  virtual real_t gm1_from_rho(real_t rho) const =0;

  /**
  @param gm1 Pseudo enthalpy \f$ g-1 \f$
  @return Rest mass density \f$ \rho \f$
  **/
  virtual real_t rho(real_t gm1) const =0;


  /**
  @param gm1 Pseudo enthalpy \f$ g-1 \f$
  @return Specific internal energy \f$\epsilon \f$
  **/
  virtual real_t eps(real_t gm1) const =0; 


  /**
  @param gm1 Pseudo enthalpy \f$ g-1 \f$
  @return Pressure \f$ P \f$ 
  **/
  virtual real_t press(real_t gm1) const =0;


  /**
  @param gm1 Pseudo enthalpy \f$ g-1 \f$
  @return Specific enthalpy \f$ h-1  \f$
  **/
  virtual real_t hm1(real_t gm1) const =0;

  /**
  @param gm1 Pseudo enthalpy \f$ g-1 \f$
  @return Adiabatic soundspeed \f$ c_s \f$
  **/
  virtual real_t csnd(real_t gm1) const =0;

  /**
  @param gm1 Pseudo enthalpy \f$ g-1 \f$
  @return Temperature \f$ T \f$ 
  
  @throws std::runtime_error if temperature is not implemented
  **/
  virtual real_t temp(real_t gm1) const =0;

  /**
  @param gm1 Pseudo enthalpy \f$ g-1 \f$
  @return Electron fraction \f$ Y_e \f$ 
  @throws std::runtime_error if electron fraction is not implemented
  **/
  virtual real_t ye(real_t gm1) const =0;

};

}// namespace implementations
}// namespace EOS_Toolkit


#endif

