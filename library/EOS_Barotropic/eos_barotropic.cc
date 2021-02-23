#include "eos_barotropic.h"
#include <stdexcept>
#include <cassert>
#include <limits>

using namespace std;
using namespace EOS_Toolkit;
using namespace EOS_Toolkit::detail;
using namespace EOS_Toolkit::implementations;


namespace {
class eos_barotr_invalid: public eos_barotr_impl {
  public:

  static std::runtime_error nope() 
  {
    return std::runtime_error("eos_barotr: uninitialized use");
  }
  
  static std::runtime_error invalid() 
  {
    return std::runtime_error(
            "eos_barotr called on invalid matter state");
  }


  eos_barotr_invalid() = default;
  [[ noreturn ]] const range& range_rho() const final {throw(nope());}
  [[ noreturn ]] const range& range_gm1() const final {throw(nope());}
  [[ noreturn ]] real_t minimal_h() const final {throw(nope());}
  [[ noreturn ]] bool is_isentropic() const final  {throw(nope());}
  [[ noreturn ]] bool is_zero_temp() const final  {throw(nope());}
  [[ noreturn ]] bool has_temp() const final  {throw(nope());}
  [[ noreturn ]] bool has_efrac() const final  {throw(nope());}

  [[ noreturn ]] 
  real_t gm1_from_rho(real_t) const final {throw(invalid());}
  [[ noreturn ]] real_t rho(real_t) const final {throw(invalid());}
  [[ noreturn ]] real_t eps(real_t) const  final {throw(invalid());}
  [[ noreturn ]] real_t press(real_t) const final {throw(invalid());}
  [[ noreturn ]] real_t hm1(real_t) const final {throw(invalid());}
  [[ noreturn ]] real_t csnd(real_t) const final {throw(invalid());}
  [[ noreturn ]] real_t temp(real_t) const final {throw(invalid());}
  [[ noreturn ]] real_t ye(real_t) const final {throw(invalid());}
};

}


const eos_barotr_base::spimpl_t eos_barotr_base::pbad{ 
  std::make_shared<eos_barotr_invalid>()
};



auto eos_barotr::is_rho_valid(real_t rho) const -> bool
{
  return impl().range_rho().contains(rho);
}

auto eos_barotr::is_gm1_valid(real_t gm1) const -> bool
{
  return impl().range_gm1().contains(gm1);
}


auto eos_barotr::at_rho(real_t rho) const -> state
{
  if (!is_rho_valid(rho)) return {};
  return {impl(), impl().gm1_from_rho(rho)};
}

auto eos_barotr::at_gm1(real_t gm1) const -> state
{
  if (!is_gm1_valid(gm1)) return {};
  return {impl(), gm1};
}



auto eos_barotr::state::gm1() const -> real_t
{ 
  if (!am_ok()) throw(eos_barotr_invalid::invalid());
  assert(gm1_ >= 0);
  return gm1_;
}


auto eos_barotr::state::press() const -> real_t 
{
  real_t press = impl().press(gm1_);  
  assert(press >= 0);
  return press;
}

auto eos_barotr::state::rho() const -> real_t 
{
  real_t rho = impl().rho(gm1_);  
  assert(rho >= 0);
  return rho;
}

auto eos_barotr::state::csnd() const -> real_t 
{
  real_t cs = impl().csnd(gm1_);  
  assert(cs < 1.0);
  assert(cs >= 0);
  return cs;
}

auto eos_barotr::state::temp() const -> real_t 
{
  real_t temp = impl().temp(gm1_);  
  assert(temp >= 0);
  return temp;
}

auto eos_barotr::state::eps() const -> real_t
{
  real_t eps = impl().eps(gm1_);  
  assert(eps >= -1);
  return eps;
}

auto eos_barotr::state::hm1() const -> real_t
{
  real_t hm1 = impl().hm1(gm1_);  
  assert(hm1 > -1);
  return hm1;
}


auto eos_barotr::state::ye() const -> real_t
{
  return impl().ye(gm1_);
}

auto eos_barotr::gm1_at_rho(real_t rho) const -> real_t
{
  auto s = at_rho(rho);
  return s ? s.gm1() : numeric_limits<real_t>::quiet_NaN();
}

auto eos_barotr::press_at_rho(real_t rho) const -> real_t
{
  auto s = at_rho(rho);
  return s ? s.press() : numeric_limits<real_t>::quiet_NaN();
}

auto eos_barotr::eps_at_rho(real_t rho) const -> real_t
{
  auto s = at_rho(rho);
  return s ? s.eps() : numeric_limits<real_t>::quiet_NaN();
}

auto eos_barotr::hm1_at_rho(real_t rho) const -> real_t
{
  auto s = at_rho(rho);
  return s ? s.hm1() : numeric_limits<real_t>::quiet_NaN();
}

auto eos_barotr::csnd_at_rho(real_t rho) const -> real_t
{
  auto s = at_rho(rho);
  return s ? s.csnd() : numeric_limits<real_t>::quiet_NaN();
}

auto eos_barotr::temp_at_rho(real_t rho) const -> real_t
{
  auto s = at_rho(rho);
  return s ? s.temp() : numeric_limits<real_t>::quiet_NaN();
}

auto eos_barotr::ye_at_rho(real_t rho) const -> real_t
{
  auto s = at_rho(rho);
  return s ? s.ye() : numeric_limits<real_t>::quiet_NaN();
}

auto eos_barotr::rho_at_gm1(real_t gm1) const -> real_t
{
  auto s = at_gm1(gm1);
  return s ? s.rho() : numeric_limits<real_t>::quiet_NaN();
}

auto eos_barotr::press_at_gm1(real_t gm1) const -> real_t
{
  auto s = at_gm1(gm1);
  return s ? s.press() : numeric_limits<real_t>::quiet_NaN();
}

auto eos_barotr::eps_at_gm1(real_t gm1) const -> real_t
{
  auto s = at_gm1(gm1);
  return s ? s.eps() : numeric_limits<real_t>::quiet_NaN();
}

auto eos_barotr::hm1_at_gm1(real_t gm1) const -> real_t
{
  auto s = at_gm1(gm1);
  return s ? s.hm1() : numeric_limits<real_t>::quiet_NaN();
}

auto eos_barotr::csnd_at_gm1(real_t gm1) const -> real_t
{
  auto s = at_gm1(gm1);
  return s ? s.csnd() : numeric_limits<real_t>::quiet_NaN();
}

auto eos_barotr::temp_at_gm1(real_t gm1) const -> real_t
{
  auto s = at_gm1(gm1);
  return s ? s.temp() : numeric_limits<real_t>::quiet_NaN();
}

auto eos_barotr::ye_at_gm1(real_t gm1) const -> real_t
{
  auto s = at_gm1(gm1);
  return s ? s.ye() : numeric_limits<real_t>::quiet_NaN();
}

