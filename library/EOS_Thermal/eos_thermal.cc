/** \file eos_thermal.cc
\brief Implementation of thermal EOS interface
**/

#include "eos_thermal_impl.h"
#include <limits>
#include <stdexcept>
#include <cassert>

using namespace std;
using namespace EOS_Toolkit;
using namespace EOS_Toolkit::detail;
using namespace EOS_Toolkit::implementations;


namespace {

class eos_thermal_invalid : public eos_thermal_impl {
  public:
    
  static std::runtime_error nope() 
  {
    return std::runtime_error("eos_thermal: uninitialized use");
  }
  
  static std::runtime_error invalid() 
  {
    return std::runtime_error(
            "eos_thermal called on invalid matter state");
  }
    

  [[ noreturn ]] real_t therm_from_rho_eps_ye(real_t rho, 
    real_t eps, real_t ye) const final {throw invalid();}

  [[ noreturn ]] real_t therm_from_rho_temp_ye(real_t rho,
    real_t temp, real_t ye) const final {throw invalid();}

  [[ noreturn ]] real_t eps(real_t rho, 
    real_t th, real_t ye) const final {throw invalid();}
  
  [[ noreturn ]] real_t temp(real_t rho, 
    real_t th, real_t ye) const final {throw invalid();}

  [[ noreturn ]] real_t press(real_t rho, 
    real_t th, real_t ye) const final {throw invalid();}

  [[ noreturn ]] real_t csnd(real_t rho, 
    real_t th, real_t ye) const final {throw invalid();}
  
  [[ noreturn ]] real_t sentr(real_t rho, 
    real_t th, real_t ye) const final {throw invalid();}

  [[ noreturn ]] real_t dpress_drho(real_t rho, 
    real_t th, real_t ye) const final {throw invalid();}

  [[ noreturn ]] real_t dpress_deps(real_t rho, 
    real_t the, real_t ye) const final {throw invalid();}
  
  
  [[ noreturn ]] const range& range_rho() const final {throw nope();}
  [[ noreturn ]] const range& range_ye() const final {throw nope();}

  
  [[ noreturn ]] range 
  range_eps(real_t rho, real_t ye) const final {throw nope();}

  [[ noreturn ]] range 
  range_temp(real_t rho, real_t ye) const final {throw nope();}

  [[ noreturn ]] real_t minimal_h() const final {throw nope();}

};

}


const eos_thermal_base::spimpl_t eos_thermal_base::pbad { 
  std::make_shared<eos_thermal_invalid>()
};


auto eos_thermal_base::state_base::rho() const -> real_t 
{
  if (!valid()) throw eos_thermal_invalid::invalid();
  return rho_;
}

auto eos_thermal_base::state_base::ye() const -> real_t 
{
  if (!valid()) throw eos_thermal_invalid::invalid();
  return ye_;
}


auto eos_thermal::state::press() const -> real_t 
{
  real_t p = eos().press(rho(), therm(), ye());  
  assert(p >= 0);
  return p;
}

auto eos_thermal::state::csnd() const -> real_t 
{
  real_t cs = eos().csnd(rho(), therm(), ye());  
  assert(cs < 1.0);
  assert(cs >= 0);
  return cs;
}

auto eos_thermal::state::temp() const -> real_t 
{
  real_t temp = eos().temp(rho(), therm(), ye());  
  assert(temp >= 0);
  return temp;
}

auto eos_thermal::state::eps() const -> real_t
{
  real_t eps = eos().eps(rho(), therm(), ye());  
  assert(eps >= -1);
  return eps;
}

auto eos_thermal::state::sentr() const -> real_t
{
  return eos().sentr(rho(), therm(), ye());  
}

auto eos_thermal::state::dpress_drho() const -> real_t
{
  return eos().dpress_drho(rho(), therm(), ye());  
}

auto eos_thermal::state::dpress_deps() const -> real_t
{
  return eos().dpress_deps(rho(), therm(), ye());  
}


auto eos_thermal::at_rho_eps_ye(real_t rho, real_t eps, 
                             real_t ye) const
-> state
{
  if (!is_rho_eps_ye_valid(rho,eps,ye)) return {};
  return {impl(), rho, impl().therm_from_rho_eps_ye(rho,eps,ye), ye};
}

auto eos_thermal::at_rho_temp_ye(real_t rho, real_t temp, 
                              real_t ye) const
-> state
{
  if (!is_rho_temp_ye_valid(rho,temp,ye)) return {};
  return {impl(), rho, impl().therm_from_rho_temp_ye(rho,temp,ye), ye};
}


auto eos_thermal::range_rho() const -> const range&
{
  return impl().range_rho();
}

auto eos_thermal::range_ye() const -> const range& 
{
  return impl().range_ye();
}


auto eos_thermal::range_eps(real_t rho, real_t ye) const -> range 
{
  if (!is_rho_valid(rho))
    throw range_error("eos_thermal: specific energy range for "
                      "invalid density requested");
  if (!is_ye_valid(ye))
    throw range_error("eos_thermal: specific energy range for "
                      "invalid electron fraction requested");

  return impl().range_eps(rho, ye);
}

 
auto eos_thermal::range_temp(real_t rho, real_t ye) const -> range
{
  if (!is_rho_valid(rho))
    throw range_error("eos_thermal: temperature range for "
                      "invalid density requested");
  if (!is_ye_valid(ye))
    throw range_error("eos_thermal: temperature range for "
                      "invalid electron fraction requested");
  return impl().range_temp(rho, ye);
}

auto eos_thermal::minimal_h() const -> real_t
{
  real_t h0 = impl().minimal_h();
  assert(h0>0);
  return h0;
}


auto eos_thermal::is_rho_valid(real_t rho) const -> bool
{
  return range_rho().contains(rho);
}

auto eos_thermal::is_ye_valid(real_t ye) const -> bool
{
  return range_ye().contains(ye);
}

auto eos_thermal::is_rho_ye_valid(real_t rho, real_t ye) const 
-> bool
{
  return impl().range_rho().contains(rho) 
         && impl().range_ye().contains(ye);
}

auto eos_thermal::is_rho_eps_ye_valid(real_t rho, 
                            real_t eps, real_t ye) const -> bool
{
  return is_rho_ye_valid(rho, ye) 
         && impl().range_eps(rho, ye).contains(eps);
}

auto eos_thermal::is_rho_temp_ye_valid(real_t rho, 
                            real_t temp, real_t ye) const ->bool
{
  return is_rho_ye_valid(rho, ye) 
         && impl().range_temp(rho, ye).contains(temp);

}

auto eos_thermal::press_at_rho_eps_ye(real_t rho, real_t eps, 
                                      real_t ye) const -> real_t
{
  auto s = at_rho_eps_ye(rho, eps, ye);
  return s ? s.press() : numeric_limits<real_t>::quiet_NaN();  
}


auto eos_thermal::csnd_at_rho_eps_ye(real_t rho, real_t eps, 
                                      real_t ye) const -> real_t
{
  auto s = at_rho_eps_ye(rho, eps, ye);
  return s ? s.csnd() : numeric_limits<real_t>::quiet_NaN();  
}

auto eos_thermal::temp_at_rho_eps_ye(real_t rho, real_t eps, 
                                      real_t ye) const -> real_t
{
  auto s = at_rho_eps_ye(rho, eps, ye);
  return s ? s.temp() : numeric_limits<real_t>::quiet_NaN();  
}

auto eos_thermal::sentr_at_rho_eps_ye(real_t rho, real_t eps, 
                                      real_t ye) const -> real_t
{
  auto s = at_rho_eps_ye(rho, eps, ye);
  return s ? s.sentr() : numeric_limits<real_t>::quiet_NaN();  
}

auto eos_thermal::dpress_drho_at_rho_eps_ye(real_t rho, real_t eps, 
                                            real_t ye) const -> real_t
{
  auto s = at_rho_eps_ye(rho, eps, ye);
  return s ? s.dpress_drho() : numeric_limits<real_t>::quiet_NaN();  
}

auto eos_thermal::dpress_deps_at_rho_eps_ye(real_t rho, real_t eps, 
                                            real_t ye) const -> real_t
{
  auto s = at_rho_eps_ye(rho, eps, ye);
  return s ? s.dpress_deps() : numeric_limits<real_t>::quiet_NaN();  
}


auto eos_thermal::press_at_rho_temp_ye(real_t rho, real_t temp, 
                                      real_t ye) const -> real_t
{
  auto s = at_rho_temp_ye(rho, temp, ye);
  return s ? s.press() : numeric_limits<real_t>::quiet_NaN();  
}


auto eos_thermal::csnd_at_rho_temp_ye(real_t rho, real_t temp, 
                                      real_t ye) const -> real_t
{
  auto s = at_rho_temp_ye(rho, temp, ye);
  return s ? s.csnd() : numeric_limits<real_t>::quiet_NaN();  
}

auto eos_thermal::eps_at_rho_temp_ye(real_t rho, real_t temp, 
                                      real_t ye) const -> real_t
{
  auto s = at_rho_temp_ye(rho, temp, ye);
  return s ? s.eps() : numeric_limits<real_t>::quiet_NaN();  
}

auto eos_thermal::sentr_at_rho_temp_ye(real_t rho, real_t temp, 
                                      real_t ye) const -> real_t
{
  auto s = at_rho_temp_ye(rho, temp, ye);
  return s ? s.sentr() : numeric_limits<real_t>::quiet_NaN();  
}

auto eos_thermal::dpress_drho_at_rho_temp_ye(real_t rho, real_t temp, 
                                            real_t ye) const -> real_t
{
  auto s = at_rho_temp_ye(rho, temp, ye);
  return s ? s.dpress_drho() : numeric_limits<real_t>::quiet_NaN();  
}

auto eos_thermal::dpress_deps_at_rho_temp_ye(real_t rho, real_t temp, 
                                            real_t ye) const -> real_t
{
  auto s = at_rho_temp_ye(rho, temp, ye);
  return s ? s.dpress_deps() : numeric_limits<real_t>::quiet_NaN();  
}


eos_thermal_impl::~eos_thermal_impl() = default;

