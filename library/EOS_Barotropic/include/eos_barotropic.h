#ifndef EOS_BAROTROPIC_H
#define EOS_BAROTROPIC_H

#include "config.h"
#include <memory>
#include "eos_barotropic_internals.h"
#include "intervals.h"

namespace EOS_Toolkit {


/// Interface for generic barotropic EOS
class eos_barotr : detail::eos_barotr_base {
  public:

  ///Synonym for \ref interval
  using range  = interval<real_t>;
  
  using eos_barotr_base::impl_t;
  
  ///Class representing the matter state for the eos_barotr interface 
  class state : state_base {
    const real_t gm1_{0};
    
    protected:
      
    state() = default;
      
    state(const impl_t& i, real_t gm1) 
    : state_base{i}, gm1_{gm1} {}
      
    
    public:

    ///Copy constructor
    state(const state&)             = default;
    
    ///Move constructor
    state(state &&)                 = default;
    
    ///Objects not meant to be assigned
    state& operator=(const state&)  = delete;
    
    ///Objects not meant to be assigned
    state& operator=(state&&)  = delete;
    
    ~state()                        = default;
    
    bool valid() const noexcept {return am_ok();}
    
    /**\brief Conversion operator to bool
    @return If state is valid.
    **/
    explicit operator bool() const noexcept {return valid();}

    /**
    @return Pseudo enthalpy \f$ g - 1 \f$ 
    
    \post Guarantees \f$ g \ge 1 \f$ 
    \throws std::runtime_error if state is invalid
    **/
    auto gm1() const -> real_t;

    /**
    @return Mass density \f$ \rho \f$ 
    
    \post Guarantees \f$ \rho \ge 0 \f$ 
    \throws std::runtime_error if state is invalid
    **/
    auto rho() const -> real_t;

    /**
    @return Pressure \f$ P \f$ 
    
    \post Guarantees \f$ P\ge 0 \f$ 
    \throws std::runtime_error if state is invalid
    **/
    auto press() const -> real_t;

    /**
    @return Specific internal energy \f$ \epsilon \f$ 
    
    \post Guarantees \f$ \epsilon \ge -1 \f$ 
    \throws std::runtime_error if state is invalid 
    **/
    auto eps() const -> real_t;

    /**
    @return Specific enthalpy \f$ h - 1 \f$ 
    
    \post Guarantees \f$ h > 0 \f$ (*not* \f$ h >= 1 \f$) 
    \throws std::runtime_error if state is invalid
    **/
    auto hm1() const -> real_t;
    
    /**
    @return Speed of sound \f$ c_s \f$
    
    \post Guarantees \f$ 0\le c_s < 1 \f$ 
    \throws std::runtime_error if state is invalid
    **/
    auto csnd() const -> real_t;

    /**
    @return Temperature \f$ T \f$
    
    \post Guarantees \f$ T\ge 0 \f$ 
    \throws std::runtime_error if state is invalid
    \throws std::runtime_error if temperature not available for EOS
    **/
    auto temp() const -> real_t;
    
    /**
    @return Electron fraction \f$ Y_e \f$
    
    \throws std::runtime_error if state is invalid
    \throws std::runtime_error if composition not available for EOS
    **/
    auto ye() const -> real_t;
    
    friend class eos_barotr;
  };

  /**\brief Constructor from pointer to implementation.
  
  This is intended for EOS implementors. Users should use the
  functions provided by a given implementation to obtain the EOS. 
  
  @param eos Shared pointer to implementation
  **/
  explicit eos_barotr(spimpl_t eos) 
  : eos_barotr_base(std::move(eos)) {}

  /**\brief Default constructor

  Creates uninitialized EOS. Any attemt to use resulting object throws 
  an exception. One can use it after copying another EOS to it.
  */
  eos_barotr()                              = default;  
  
  ///Copy constructor.
  eos_barotr(const eos_barotr &)            = default;
  
  ///Move constructor
  eos_barotr(eos_barotr &&)                 = default;
  
  ///Assignment operator.
  eos_barotr& operator=(const eos_barotr&)  = default;
  
  ///Move assignment operator
  eos_barotr& operator=(eos_barotr&&)       = default;
  
  /**\brief Destructor
  
  The EOS implementation will be destructed when no other copy of the
  EOS is using it ("last person turns off the light").
  **/
  ~eos_barotr()                             = default;

  /**\brief Specify a matter state based on mass density
      
  @param rho  Mass density \f$ \rho \f$
  @returns    Object representing matter \ref state 
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto at_rho(real_t rho) const -> state;

  /**\brief Specify a matter state based on pseudo enthalpy
      
  @param gm1  Pseudo enthalpy \f$ g - 1 \f$
  @returns    Object representing matter \ref state 
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto at_gm1(real_t gm1) const -> state;

  /**
  @return Whether EOS is isentropic
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto is_isentropic() const -> bool {
    return impl().is_isentropic();
  }

  /**
  @return Whether EOS is zero-temperature (or false if temperature 
          not implemented)
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto is_zero_temp() const -> bool {
    return impl().is_zero_temp();
  }

  /**
  @return Whether EOS provides temperature
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto has_temp() const -> bool {
    return impl().has_temp();
  }

  /**
  @return Whether EOS provides electron fraction
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto has_efrac() const -> bool {
    return impl().has_efrac();
  }
  
  /**
  @returns Validity \ref range for mass density
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto range_rho() const -> const range& {
    return impl().range_rho();
  }
  
  /**
  @returns Validity \ref range for pseudo enthalpy \f$ g - 1 \f$
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto range_gm1() const -> const range& {
    return impl().range_gm1();
  }
  
  /**
  @return  Global lower bound for relativistic enthalpy.
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto minimal_h() const -> real_t {
    return impl().minimal_h();
  }
  
  /**
  @param rho  Mass density \f$ \rho \f$
  @return Wether density is in valid range
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto is_rho_valid(real_t rho) const -> bool;
  
  /**
  @param gm1 Pseudo enthalpy \f$ g - 1 \f$
  @return Whether pseudo enthalpy is in valid range
  
  \throws std::runtime_error if called for unitialized object
  **/
  auto is_gm1_valid(real_t gm1) const -> bool;


  /**\brief Compute pseudo enthalpy from mass density
      
  @param rho  Mass density \f$ \rho \f$
  @return Pseudo enthalpy \f$ g - 1 \f$ if state is valid else NAN

  \post Guarantees \f$ g \ge 1 \f$  (unless NAN)
  \throws std::runtime_error if called for unitialized object
  **/
  auto gm1_at_rho(real_t rho) const -> real_t;


  /**\brief Compute pressure from mass density
      
  @param rho  Mass density \f$ \rho \f$
  @returns    Pressure if state is valid else NAN 
  
  \post Guarantees \f$ P\ge 0 \f$  (unless NAN)
  \throws std::runtime_error if called for unitialized object
  **/
  auto press_at_rho(real_t rho) const -> real_t;
  
  /**\brief Compute specific energy from mass density
      
  @param rho  Mass density \f$ \rho \f$
  @returns    \f$ \epsilon \f$ if state is valid else NAN
  
  \post Guarantees \f$ \epsilon \ge -1 \f$  (unless NAN)
  \throws std::runtime_error if called for unitialized object
  **/
  auto eps_at_rho(real_t rho) const -> real_t;
  
  /**\brief Compute specific enthalpy from mass density
      
  @param rho  Mass density \f$ \rho \f$
  @returns    \f$ h - 1 \f$ if state is valid else NAN
  
  \post Guarantees \f$ h > 0 \f$  (unless NAN)
  \throws std::runtime_error if called for unitialized object
  **/
  auto hm1_at_rho(real_t rho) const -> real_t;
  
  /**\brief Compute sound speed from mass density
      
  @param rho  Mass density \f$ \rho \f$
  @returns    \f$ c_s \f$ if state is valid else NAN
  
  \post Guarantees \f$ 0\le c_s < 1 \f$  (unless NAN)
  \throws std::runtime_error if called for unitialized object
  **/
  auto csnd_at_rho(real_t rho) const -> real_t;
  
  /**\brief Compute temperature from mass density
      
  @param rho  Mass density \f$ \rho \f$
  @returns    \f$ T \f$ if state is valid else NAN
  
  \post Guarantees \f$ T\ge 0 \f$ (unless NAN)
  \throws std::runtime_error if temperature not available for EOS
  \throws std::runtime_error if called for unitialized object
  **/
  auto temp_at_rho(real_t rho) const -> real_t;
  
  /**\brief Compute electron fraction from mass density
      
  @param rho  Mass density \f$ \rho \f$
  @returns    \f$ Y_e \f$ if state is valid else NAN
  
  \throws std::runtime_error if electron fraction not available for EOS
  \throws std::runtime_error if called for unitialized object
  **/
  auto ye_at_rho(real_t rho) const -> real_t;
  
  
  /**\brief Compute mass density from pseudo enthalpy 
      
  @param gm1 Pseudo enthalpy \f$ g - 1 \f$ 
  @return \f$ \rho \f$ if state is valid else NAN

  \post Guarantees \f$ \rho \ge 0 \f$ (unless NAN)
  \throws std::runtime_error if called for unitialized object
  **/
  auto rho_at_gm1(real_t gm1) const -> real_t;

  /**\brief Compute pressure from pseudo enthalpy
      
  @param gm1 Pseudo enthalpy \f$ g - 1 \f$ 
  @returns    Pressure if state is valid else NAN 
  
  \post Guarantees \f$ P\ge 0 \f$ (unless NAN)
  \throws std::runtime_error if called for unitialized object
  **/
  auto press_at_gm1(real_t gm1) const -> real_t;
  
  /**\brief Compute specific energy from pseudo enthalpy
      
  @param gm1 Pseudo enthalpy \f$ g - 1 \f$ 
  @returns    \f$ \epsilon \f$ if state is valid else NAN
  
  \post Guarantees \f$ \epsilon \ge -1 \f$  (unless NAN)
  \throws std::runtime_error if called for unitialized object
  **/
  auto eps_at_gm1(real_t gm1) const -> real_t;
  
  /**\brief Compute specific enthalpy  from pseudo enthalpy
      
  @param gm1 Pseudo enthalpy \f$ g - 1 \f$ 
  @returns    \f$ h - 1 \f$ if state is valid else NAN
  
  \post Guarantees \f$ h > 0 \f$  (unless NAN)
  \throws std::runtime_error if called for unitialized object
  **/
  auto hm1_at_gm1(real_t gm1) const -> real_t;
  
  /**\brief Compute sound speed  from pseudo enthalpy
      
  @param gm1 Pseudo enthalpy \f$ g - 1 \f$ 
  @returns    \f$ c_s \f$ if state is valid else NAN
  
  \post Guarantees \f$ 0\le c_s < 1 \f$ (unless NAN)
  \throws std::runtime_error if called for unitialized object
  **/
  auto csnd_at_gm1(real_t gm1) const -> real_t;
  
  /**\brief Compute temperature from  from pseudo enthalpy
      
  @param gm1 Pseudo enthalpy \f$ g - 1 \f$ 
  @returns    \f$ T \f$ if state is valid else NAN
  
  \post Guarantees \f$ T\ge 0 \f$ (unless NAN)
  \throws std::runtime_error if temperature not available for EOS
  \throws std::runtime_error if called for unitialized object
  **/
  auto temp_at_gm1(real_t gm1) const -> real_t;
  
  /**\brief Compute electron fraction  from pseudo enthalpy
      
  @param gm1 Pseudo enthalpy \f$ g - 1 \f$ 
  @returns    \f$ Y_e \f$ if state is valid else NAN
  
  \throws std::runtime_error if electron fraction not available for EOS
  \throws std::runtime_error if called for unitialized object
  **/
  auto ye_at_gm1(real_t gm1) const -> real_t;

};





}// namespace EOS_Toolkit


#endif

