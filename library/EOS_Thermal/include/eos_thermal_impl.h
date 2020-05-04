/*! \file eos_thermal_impl.h
\brief Defines generic API for thermal EOS implementations.
Defines class eos_thermal_impl, the abstract interface from which 
all EOS implementations are derived.
*/
#ifndef EOSTHERMAL_IMPL_H
#define EOSTHERMAL_IMPL_H

#include "eos_thermal.h"


namespace EOS_Toolkit {
namespace implementations {

/// Virtual abstract base class defining Thermal EOS implementation API.
class eos_thermal_impl {
  public:
  using range  = eos_thermal::range;

  //EOS implementations are not meant to be copied
  eos_thermal_impl(const eos_thermal_impl &)            = delete;
  
  //EOS implementations are not meant to be moved
  eos_thermal_impl(eos_thermal_impl &&)                 = delete;
  
  //EOS implementations are not meant to be assigned to
  eos_thermal_impl& operator=(const eos_thermal_impl &) = delete;
  
  //EOS implementations are not meant to be move-assigned to
  eos_thermal_impl& operator=(eos_thermal_impl &&)      = delete;
  
  //Virtual destructor
  virtual ~eos_thermal_impl();

  /**
  @return Internal representation \f$ \theta \f$ of the 
          thermal degree of freedom.
          
  
  @param rho Rest mass density  \f$ \rho \f$
  @param eps Specific internal energy \f$ \epsilon \f$
  @param ye  Electron fraction \f$ Y_e \f$
  **/
  virtual real_t therm_from_rho_eps_ye(real_t rho, 
                                       real_t eps, real_t ye) const=0;

  /**
  @return Internal representation \f$ \theta \f$ of the 
          thermal degree of freedom.
  
  @param rho Rest mass density  \f$ \rho \f$
  @param temp Temperature \f$ T \f$
  @param ye  Electron fraction \f$ Y_e \f$
  **/
  virtual real_t therm_from_rho_temp_ye(real_t rho, 
                                        real_t temp, real_t ye) const=0;
  
  /**
  @return Pressure \f$ P \f$
  
  @param rho   Rest mass density  \f$ \rho \f$
  @param therm Thermal variable \f$ \theta \f$
  @param ye    Electron fraction \f$ Y_e \f$
  **/
  virtual real_t press(real_t rho, real_t therm, real_t ye) const=0;
  
  
  /**
  @return Adiabatic soundspeed \f$ c_s \f$
  
  @param rho   Rest mass density  \f$ \rho \f$
  @param therm Thermal variable \f$ \theta \f$
  @param ye    Electron fraction \f$ Y_e \f$
  **/
  virtual real_t csnd(real_t rho, real_t therm, real_t ye) const=0;
  
  /**
  @return Specific internal energy \f$ \epsilon \f$
  
  @param rho   Rest mass density  \f$ \rho \f$
  @param therm Thermal variable \f$ \theta \f$
  @param ye    Electron fraction \f$ Y_e \f$
  **/
  virtual real_t eps(real_t rho, real_t therm, real_t ye) const=0;
  
  /**
  @return Temperature \f$ T \f$
  
  @param rho   Rest mass density  \f$ \rho \f$
  @param therm Thermal variable \f$ \theta \f$
  @param ye    Electron fraction \f$ Y_e \f$
  **/
  virtual real_t temp(real_t rho, real_t therm, real_t ye) const=0;
  
  /**
  @return Specific entropy \f$ s \f$
  
  @param rho   Rest mass density  \f$ \rho \f$
  @param therm Thermal variable \f$ \theta \f$
  @param ye    Electron fraction \f$ Y_e \f$
  **/
  virtual real_t sentr(real_t rho, real_t therm, real_t ye) const=0;

  
  /**
  @return Partial derivative 
          \f$ \frac{\partial P}{\partial \rho} \f$
  
  @param rho   Rest mass density  \f$ \rho \f$
  @param therm Thermal variable \f$ \theta \f$
  @param ye    Electron fraction \f$ Y_e \f$
  **/
  virtual real_t dpress_drho(real_t rho, 
                             real_t therm, real_t ye) const=0;

  /**
  @return Partial derivative 
          \f$ \frac{\partial P}{\partial \epsilon} \f$
  
  @param rho   Rest mass density  \f$ \rho \f$
  @param therm Thermal variable \f$ \theta \f$
  @param ye    Electron fraction \f$ Y_e \f$
  **/
  virtual real_t dpress_deps(real_t rho, 
                             real_t therm, real_t ye) const=0;
  
  
  /**
  @return Valid range for density \f$ \rho \f$
  **/
  virtual auto range_rho() const -> const range& = 0;

  /**
  @return Valid range for electron fraction \f$ Y_e \f$
  **/
  virtual auto range_ye() const -> const range& = 0;

  /** 
  @return Valid range for specifc energy \f$ \epsilon \f$
  
  @param rho   Rest mass density  \f$ \rho \f$
  @param ye    Electron fraction \f$ Y_e \f$
  **/
  virtual range range_eps(real_t rho, real_t ye) const=0;

  
  /** 
  @return Valid range for temperature \f$ T \f$
  
  @param rho   Rest mass density  \f$ \rho \f$
  @param ye    Electron fraction \f$ Y_e \f$
  **/
  virtual range range_temp(real_t rho, real_t ye) const=0;


  /**
  @return Global minimum \f$ h_0 \f$ of enthalpy.
  
  It must be guaranteed that the enthaply computed by the EOS for 
  any valid state is always above this value. 
  
  \post The minimum must be strictly positive.
  **/
  virtual real_t minimal_h() const=0;


  protected:
  
  eos_thermal_impl() = default;


};


} // namespace implementations



} // namespace EOS_Toolkit 


#endif
