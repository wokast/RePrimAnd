#include "eos_idealgas.h"
#include "eos_idealgas_impl.h"
#include <stdexcept>
#include <algorithm>
#include <cmath>

using namespace EOS_Toolkit;
using namespace EOS_Toolkit::implementations;
using namespace std;

eos_idealgas::eos_idealgas(real_t n_, real_t max_eps_, real_t max_rho_) 
: gamma{1.0 + 1.0/n_}, gm1{1.0/n_}, rgrho{0., max_rho_}, 
  rgye{0., 1.}, min_h{1.}
{
  if (n_ < 0) {
    throw runtime_error("eos_idealgas: initialized with gamma < 1");
  } 
  if (gamma > 2.0) { // Ensure subluminal Soundspeed and P < E
    max_eps_ = min(max_eps_, 1.0 / (gamma*(gamma - 2.0)));
  }
  rgeps = {0., max_eps_};
}


/**
Using the formula 
\f[ P = \left(\Gamma - 1\right) \rho \epsilon \f]
**/
real_t eos_idealgas::press(real_t rho, real_t eps, real_t ye) const
{
  return gm1 * rho * eps;
}

/**
Using the formula 
\f[ c_s^2 = \frac{\left(\Gamma - 1\right) \epsilon}{\epsilon + 1/\Gamma} \f]
**/
real_t eos_idealgas::csnd(real_t rho, real_t eps, real_t ye) const
{
  return sqrt(gm1 * eps / (eps + 1.0/gamma));
}

real_t eos_idealgas::temp(real_t rho, real_t eps, real_t ye) const
{
  throw logic_error("eos_idealgas: temperature not implemented");
} 

real_t eos_idealgas::dpress_drho(real_t rho, real_t eps, 
                               real_t ye) const
{
  return gm1 * eps;
}

real_t eos_idealgas::dpress_deps(real_t rho, real_t eps, 
                               real_t ye) const
{
  return gm1 * rho;
}

real_t eos_idealgas::sentr(real_t rho, real_t eps, real_t ye) const
{
  throw logic_error("eos_idealgas: entropy not implemented");
}

real_t eos_idealgas::therm_from_rho_temp_ye(real_t rho, 
                                real_t temp, real_t ye) const
{
  throw logic_error("eos_idealgas: temperature not implemented");
}

auto eos_idealgas::range_eps(real_t rho, real_t ye) const -> range
{
  return rgeps;
}

auto eos_idealgas::range_temp(real_t rho, real_t ye) const -> range
{
  throw logic_error("eos_idealgas: temperature not implemented");
}

eos_thermal EOS_Toolkit::make_eos_idealgas(real_t n,  
    real_t max_eps, real_t max_rho)
{
  return eos_thermal{make_shared<eos_idealgas>(n, max_eps, max_rho)};
}

