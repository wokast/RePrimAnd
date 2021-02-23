#ifndef EOS_THERMAL_DETAILS_H
#define EOS_THERMAL_DETAILS_H

#include <memory>

#include "config.h"
#include "intervals.h"

namespace EOS_Toolkit {
  

class eos_thermal;


namespace implementations {

  class eos_thermal_impl;

}

namespace detail {

class eos_thermal_base {
  public:

  using impl_t = implementations::eos_thermal_impl;
  using spimpl_t = std::shared_ptr<const impl_t>;



  class state_base {
    const impl_t& p{eos_thermal_base::ibad()};
    const bool ok{false};


    const real_t rho_{0};
    const real_t therm_{0};
    const real_t ye_{0};

    protected:

    /// Obtain thermal variable
    /** This variable is intended only for internal use **/
    auto therm() const -> real_t {return therm_;}
      
    auto eos()   const -> const impl_t& {return p;}

    state_base() = default;
    state_base(const impl_t& i, real_t rho, 
                      real_t therm, real_t ye)
    : p(i), ok{true}, rho_(rho), therm_(therm), ye_(ye) {}
    

    public:

    state_base(const state_base&) = default;
    state_base(state_base &&) = default;

    state_base& operator=(const state_base&) = delete;
    state_base& operator=(state_base&&) = delete;
    ~state_base() = default;
    
    /// Whether state is valid
    auto valid() const -> bool {return ok;}
      
    /// Obtain mass density
    auto rho() const -> real_t;

    /// Obtain electron fraction
    auto ye() const -> real_t;
    
  };

    
  eos_thermal_base() = default;

  explicit eos_thermal_base(
    spimpl_t eosp   ///< Shared pointer to implementation
  ) : pimpl(std::move(eosp)) {}

  eos_thermal_base(const eos_thermal_base&)            = default;
  eos_thermal_base& operator=(const eos_thermal_base&) = default;
  eos_thermal_base(eos_thermal_base&&)                 = default;
  eos_thermal_base& operator=(eos_thermal_base&&)      = default;
  ~eos_thermal_base()                                  = default;

  private:  
  
  const static spimpl_t pbad; ///< Dummy implementation that always throws

  spimpl_t pimpl{pbad};  ///< Shared pointer to implementation
  
  protected:
  
  auto impl() const -> const impl_t& {return *pimpl;}          
  static auto ibad() -> const impl_t& {return *pbad;}          

};




} //namespace detail





}


#endif
