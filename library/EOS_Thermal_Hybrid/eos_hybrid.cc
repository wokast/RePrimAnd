#include "eos_hybrid.h"
#include "eos_hybrid_impl.h"

#include <stdexcept>
#include <algorithm>
#include <cmath>

using namespace EOS_Toolkit;
using namespace EOS_Toolkit::implementations;
using namespace std;


eos_hybrid::eos_hybrid(eos_barotr eos_c_,  real_t gamma_th_, 
                  real_t eps_max_, real_t rho_max_)
: eos_c{eos_c_}, gamma_th{gamma_th_}, gm1_th{gamma_th - 1}, 
  eps_max{eps_max_}, rgrho{0.0, rho_max_}, rgye{0.,1.}, 
  min_h{eos_c_.minimal_h()}
{
  if (!eos_c.is_zero_temp())
    throw runtime_error("eos_hybrid: cold EOS not zero temperature");
  if ((gamma_th>2) || (gamma_th<=1)) {
    throw runtime_error("eos_hybrid: gamma_thermal must be in the "
                        "range (1,2]");
  }  
}


real_t eos_hybrid::eps_cold(real_t rho) const
{
  return eos_c.at_rho(rho).eps();
}

real_t eos_hybrid::p_cold(real_t rho) const
{
  return eos_c.at_rho(rho).press();
}

real_t eos_hybrid::hm1_cold(real_t rho) const
{
  return eos_c.at_rho(rho).hm1();
}

real_t eos_hybrid::cs2_cold(real_t rho) const
{
  real_t cs = eos_c.at_rho(rho).csnd();
  return cs*cs;
}



real_t eos_hybrid::press(real_t rho, real_t eps, real_t ye) const
{
  real_t p_c   = p_cold(rho);
  real_t eps_c = eps_cold(rho);
  real_t p_th  = gm1_th * rho * (eps - eps_c);
  return p_c + p_th;
}

real_t eos_hybrid::csnd(real_t rho, real_t eps, real_t ye) const
{
  real_t cs2_c   = cs2_cold(rho);
  real_t eps_c   = eps_cold(rho);
  real_t h_c     = 1.0 + hm1_cold(rho);
  real_t eps_th  = eps - eps_c;
  real_t h_th    = gamma_th * eps_th;
  real_t w       = h_th / (h_c + h_th);
  real_t cs2     = (1.0 - w) * cs2_c + w * gm1_th;
  return sqrt(cs2);
}

real_t eos_hybrid::temp(real_t rho, real_t eps, real_t ye) const
{
  throw runtime_error("eos_hybrid: temperature not implemented");
} 

real_t eos_hybrid::dpress_drho(real_t rho, real_t eps, real_t ye) const
{
  real_t p_c     = p_cold(rho);
  real_t cs2_c   = cs2_cold(rho);
  real_t eps_c   = eps_cold(rho);
  real_t h_c     = 1.0 + hm1_cold(rho);
  real_t eps_th  = eps - eps_c;
  return h_c * cs2_c + gm1_th * (eps_th - p_c / rho);
}

real_t eos_hybrid::dpress_deps(real_t rho, real_t eps, real_t ye) const
{
  return gm1_th * rho;
}

real_t eos_hybrid::sentr(real_t rho, real_t eps, real_t ye) const
{
  throw runtime_error("eos_hybrid: entropy not implemented");
}

real_t eos_hybrid::therm_from_rho_temp_ye(real_t rho, 
            real_t temp, real_t ye) const
{
  throw runtime_error("eos_hybrid: temperature not implemented");
}

eos_thermal_impl::range 
eos_hybrid::range_eps(real_t rho, real_t ye) const
{
  real_t eps_min = eps_cold(rho);
  return {eps_min, eps_max};
}

eos_thermal_impl::range 
eos_hybrid::range_temp(real_t rho, real_t ye) const
{
  throw runtime_error("eos_hybrid: temperature not implemented");
}

eos_thermal EOS_Toolkit::make_eos_hybrid(EOS_Toolkit::eos_barotr eos_c,  
             real_t gamma_th, real_t eps_max, real_t rho_max)
{
  return eos_thermal{std::make_shared<eos_hybrid>(eos_c, gamma_th, 
                                                  eps_max, rho_max)};
}



