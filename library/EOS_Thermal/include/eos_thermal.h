/*! \file eos_thermal.h
\brief Defines a generic interface for thermal EOS
*/

#ifndef EOSTHERMAL_H
#define EOSTHERMAL_H

#include "eos_thermal_internals.h"

#include <string>



namespace EOS_Toolkit {


/// Interface for generic thermal EOS
class eos_thermal : detail::eos_thermal_base {

  public:
  
  ///Synonym for \ref interval
  using range  = interval<real_t>;
  
  

  ///Class representing the matter state for the eos_thermal interface 
  class state : public state_base {
    public:
    
    using state_base::state_base;
    
    /**\brief Conversion operator to bool
    @return If state is valid.
    **/
    explicit operator bool() const {return valid();}
    
    /**
    @return Pressure \f$ P \f$ 
    
    \post Guarantees \f$ P\ge 0 \f$ 
    \throws std::runtime_error if state is invalid
    **/
    auto press() const -> real_t;
    
    /**
    @return Speed of sound \f$ c_s \f$
    
    \post Guarantees \f$ 0\le c_s < 1 \f$ 
    \throws std::runtime_error if state is invalid
    **/
    auto csnd() const -> real_t;

    /**
    @return Temperature \f$ T\ge 0 \f$
    
    \post Guarantees \f$ T\ge 0 \f$ 
    \throws std::runtime_error if state is invalid
    \throws std::runtime_error if temperature not available for EOS
    **/
    auto temp() const -> real_t;
    
    /**
    @return Specific entropy \f$ s \f$ 
    
    \throws std::runtime_error if state is invalid
    \throws std::runtime_error if entropy not available for EOS
    **/
    auto sentr() const -> real_t;
    
    /**
    @return Specific internal energy \f$ \epsilon \f$ 
    
    \post Guarantees \f$ \epsilon \ge -1 \f$ 
    \throws std::runtime_error if state is invalid 
    **/
    auto eps() const -> real_t;
    
    /**
    @return partial derivative of pressure with respect to 
    mass density 
    \f$ \frac{\partial P}{\partial \rho} \f$
    
    \throws std::runtime_error if state is invalid 
    **/
    auto dpress_drho() const -> real_t;
    
    /**
    @return  partial derivative of pressure with respect to 
    specific energy
    \f$ \frac{\partial P}{\partial \epsilon} \f$
        
    \throws std::runtime_error if state is invalid 
    **/
    auto dpress_deps() const -> real_t;
    
    friend class eos_thermal;
  };

    
  
    
  /**\brief Default constructor

  Creates uninitialized EOS. Any attemt to use resulting object throws 
  an exception. One can use it after copying another EOS to it.
  */
  eos_thermal() = default;

  /**\brief Constructor from pointer to implementation.
  
  This is intended for EOS implementors. Users should use the
  functions provided by a given implementation to obtain the EOS. 
  @param eosp Shared pointer to implementation
  **/
  explicit eos_thermal(spimpl_t eosp) 
  : eos_thermal_base(std::move(eosp)) {}

  ///Copy constructor.
  eos_thermal(const eos_thermal&) = default;

  ///Assignment operator.
  eos_thermal& operator=(const eos_thermal&) = default;

  ///Move constructor
  eos_thermal(eos_thermal&&) = default;
  
  ///Move assignment operator
  eos_thermal& operator=(eos_thermal&&) = default;
  
  /**\brief Destructor
  
  The EOS implementation will be destructed when no other copy of the
  EOS is using it ("last person turns off the light").
  **/
  ~eos_thermal() = default;

  /**\brief Specify a matter state based on density, specific energy, 
  and electron fraction
      
  @param rho  Mass density \f$ \rho \f$
  @param eps  Specific internal energy \f$ \epsilon \f$
  @param ye   Electron fraction \f$ Y_e \f$
  @returns    Object representing matter \ref state 
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto at_rho_eps_ye(real_t rho, real_t eps, real_t ye) const 
  -> state;

  /**\brief Specify a matter state based on density, temperature, 
  and electron fraction

  @param rho  Mass density \f$ \rho \f$
  @param temp Temperature \f$ T \f$
  @param ye   Electron fraction \f$ Y_e \f$
  @returns    Object representing matter \ref state 
  
  \throws std::runtime_error if EOS does not support temperature
  \throws std::runtime_error if called for unitialized object
  **/
  auto at_rho_temp_ye(real_t rho, real_t temp, real_t ye) const 
  -> state;
                          
  /** 
  @return validity \ref range for density
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto range_rho() const -> const range&;

  /**
  @returns validity \ref range for electron fraction
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto range_ye() const -> const range&;

  /**
  @param rho  Mass density \f$ \rho \f$
  @param ye   Electron fraction \f$ Y_e \f$
  @return validity range for specific energy.
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto range_eps(real_t rho, real_t ye) const -> range;

  /**
  @param rho  Mass density \f$ \rho \f$
  @param ye   Electron fraction \f$ Y_e \f$
  @return validity range for temperature
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto range_temp(real_t rho, real_t ye) const -> range; 
  
  /**
  @return  global lower bound for relativistic enthalpy.
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto minimal_h() const -> real_t;

  /**
  @param rho  Mass density \f$ \rho \f$
  @return if density is in valid range
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto is_rho_valid(real_t rho) const -> bool;
  
  /**
  @param ye   Electron fraction \f$ Y_e \f$
  @return if electron fraction is in valid range
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto is_ye_valid(real_t ye) const -> bool;

  /**
  @param rho  Mass density \f$ \rho \f$
  @param ye   Electron fraction \f$ Y_e \f$
  @return if density and electron fraction are in valid range
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto is_rho_ye_valid(real_t rho, real_t ye) const -> bool;

  /**
  @param rho  Mass density \f$ \rho \f$
  @param eps  Specific internal energy \f$ \epsilon \f$
  @param ye   Electron fraction \f$ Y_e \f$
  @return if density, specfic energy, and electron fraction are valid.
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto is_rho_eps_ye_valid(real_t rho, 
                           real_t eps, real_t ye) const -> bool;

  /**
  @param rho  Mass density \f$ \rho \f$
  @param temp temperature \f$ T \f$
  @param ye   Electron fraction \f$ Y_e \f$
  @return if density, temperature, and electron fraction are valid.
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto is_rho_temp_ye_valid(real_t rho, 
                            real_t temp, real_t ye) const -> bool;

  /**\brief Compute pressure based on density, specific energy, 
  and electron fraction
      
  @param rho  Mass density \f$ \rho \f$
  @param eps  Specific internal energy \f$ \epsilon \f$
  @param ye   Electron fraction \f$ Y_e \f$
  @returns    Pressure \f$ P \f$ if state is valid, else NAN 
  
  \post Guarantees \f$ P\ge 0 \f$ 
    \throws std::runtime_error if called for unitialized object
  **/
  auto press_at_rho_eps_ye(real_t rho, real_t eps, real_t ye) const 
  -> real_t;

  /**\brief Compute soundspeed based on density, specific energy, 
  and electron fraction
      
  @param rho  Mass density \f$ \rho \f$
  @param eps  Specific internal energy \f$ \epsilon \f$
  @param ye   Electron fraction \f$ Y_e \f$
  @return Speed of sound \f$ c_s \f$ if state is valid, else NAN

  \post Guarantees \f$ 0\le c_s < 1 \f$ 
  \throws std::runtime_error if called for unitialized object
  **/
  auto csnd_at_rho_eps_ye(real_t rho, real_t eps, real_t ye) const 
  -> real_t;

  /**\brief Compute temperature based on density, specific energy, 
  and electron fraction
      
  @param rho  Mass density \f$ \rho \f$
  @param eps  Specific internal energy \f$ \epsilon \f$
  @param ye   Electron fraction \f$ Y_e \f$
  @returns    Temperature \f$ T \f$ if state is valid, else NAN 
  
  \post Guarantees \f$ T\ge 0 \f$ 
  \throws std::runtime_error if temperature not available for EOS
  \throws std::runtime_error if called for unitialized object
  **/
  auto temp_at_rho_eps_ye(real_t rho, real_t eps, real_t ye) const 
  -> real_t;

  /**\brief Compute specific entropy based on density, specific energy, 
  and electron fraction
      
  @param rho  Mass density \f$ \rho \f$
  @param eps  Specific internal energy \f$ \epsilon \f$
  @param ye   Electron fraction \f$ Y_e \f$
  @returns    Specific entropy \f$ s \f$ if state is valid, else NAN 
    
  \throws std::runtime_error if entropy not available for EOS
  \throws std::runtime_error if called for unitialized object
  **/
  auto sentr_at_rho_eps_ye(real_t rho, real_t eps, real_t ye) const 
  -> real_t;

  /**\brief Compute partial derivative of pressure with respect to 
  mass density, \f$ \frac{\partial P}{\partial \rho} \f$,
  based on density, specific energy, and electron fraction.
      
  @param rho  Mass density \f$ \rho \f$
  @param eps  Specific internal energy \f$ \epsilon \f$
  @param ye   Electron fraction \f$ Y_e \f$
  @returns    Pressure derivative if state is valid, else NAN 
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto dpress_drho_at_rho_eps_ye(real_t rho, real_t eps, 
                                 real_t ye) const -> real_t;

  /**\brief Compute partial derivative of pressure with respect to 
  specific energy, \f$ \frac{\partial P}{\partial \epsilon} \f$,
  based on density, specific energy, and electron fraction.
      
  @param rho  Mass density \f$ \rho \f$
  @param eps  Specific internal energy \f$ \epsilon \f$
  @param ye   Electron fraction \f$ Y_e \f$
  @returns    Pressure derivative if state is valid, else NAN 
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto dpress_deps_at_rho_eps_ye(real_t rho, real_t eps, 
                                 real_t ye) const -> real_t;


  /**\brief Compute pressure based on density, temperature, 
  and electron fraction
      
  @param rho  Mass density \f$ \rho \f$
  @param temp Temperature \f$ T \f$
  @param ye   Electron fraction \f$ Y_e \f$
  @returns    Pressure \f$ P \f$ if state is valid, else NAN 
  
  \post Guarantees \f$ P\ge 0 \f$ 
  \throws std::runtime_error if temperature not available for EOS
  \throws std::runtime_error if called for unitialized object
  **/
  auto press_at_rho_temp_ye(real_t rho, real_t temp, real_t ye) const 
  -> real_t;

  /**\brief Compute soundspeed based on density, temperature, 
  and electron fraction
      
  @param rho  Mass density \f$ \rho \f$
  @param temp Temperature \f$ T \f$
  @param ye   Electron fraction \f$ Y_e \f$
  @return Speed of sound \f$ c_s \f$ if state is valid, else NAN

  \post Guarantees \f$ 0\le c_s < 1 \f$ 
  \throws std::runtime_error if temperature not available for EOS
  \throws std::runtime_error if called for unitialized object
  **/
  auto csnd_at_rho_temp_ye(real_t rho, real_t temp, real_t ye) const 
  -> real_t;

  /**\brief Compute specific energy based on density, temperature, 
  and electron fraction
      
  @param rho  Mass density \f$ \rho \f$
  @param temp Temperature \f$ T \f$
  @param ye   Electron fraction \f$ Y_e \f$
  @returns    \f$ \epsilon \f$ if state is valid, else NAN 
  
  \post Guarantees \f$ \epsilon \ge -1 \f$ 
  \throws std::runtime_error if temperature not available for EOS
  \throws std::runtime_error if called for unitialized object
  **/
  auto eps_at_rho_temp_ye(real_t rho, real_t temp, real_t ye) const 
  -> real_t;

  /**\brief Compute specific entropy based on density, temperature, 
  and electron fraction
      
  @param rho  Mass density \f$ \rho \f$
  @param temp Temperature \f$ T \f$
  @param ye   Electron fraction \f$ Y_e \f$
  @returns    Specific entropy \f$ s \f$ if state is valid, else NAN 

  \throws std::runtime_error if temperature not available for EOS    
  \throws std::runtime_error if entropy not available for EOS
  \throws std::runtime_error if called for unitialized object
  **/
  auto sentr_at_rho_temp_ye(real_t rho, real_t temp, real_t ye) const 
  -> real_t;

  /**\brief Compute partial derivative of pressure with respect to 
  mass density, \f$ \frac{\partial P}{\partial \rho} \f$,
  based on density, temperature, and electron fraction.
      
  @param rho  Mass density \f$ \rho \f$
  @param temp Temperature \f$ T \f$
  @param ye   Electron fraction \f$ Y_e \f$
  @returns    Pressure derivative if state is valid, else NAN 
  
  \throws std::runtime_error if temperature not available for EOS      
  \throws std::runtime_error if called for unitialized object
  **/
  auto dpress_drho_at_rho_temp_ye(real_t rho, real_t temp, 
                                  real_t ye) const -> real_t;

  /**\brief Compute partial derivative of pressure with respect to 
  specific energy, \f$ \frac{\partial P}{\partial \epsilon} \f$,
  based on density, temperature, and electron fraction.
      
  @param rho  Mass density \f$ \rho \f$
  @param temp Temperature \f$ T \f$
  @param ye   Electron fraction \f$ Y_e \f$
  @returns    Pressure derivative if state is valid, else NAN 
  
  \throws std::runtime_error if temperature not available for EOS      
  \throws std::runtime_error if called for unitialized object
  **/
  auto dpress_deps_at_rho_temp_ye(real_t rho, real_t temp, 
                                  real_t ye) const -> real_t;

 
};


} // namespace EOS_Toolkit





#endif

