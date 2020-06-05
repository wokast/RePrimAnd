#include "eos_barotr_gpoly.h"
#include "eos_barotr_gpoly_impl.h"
#include <cmath>
#include <stdexcept>

using namespace std;
using namespace EOS_Toolkit;
using namespace EOS_Toolkit::implementations;



eos_barotr_gpoly::eos_barotr_gpoly(real_t n_, real_t rmd_p_, 
                                   real_t sed0_, real_t rho_max_) 
: rgrho(0, rho_max_), n(n_), rmd_p(rmd_p_), np1(n + 1), 
  gamma(1.0 + 1.0 / n), invn(1.0 / n), sed0(sed0_), h0(1.0 + sed0_)
{
  rggm1 = {0, gm1_from_rho(rho_max_)};
}


/**
\return \f$ g-1 =  \frac{(n+1)}{h_0} \left(\frac{\rho}{\rho_p}\right)^{1/n} \f$
*/
real_t eos_barotr_gpoly::gm1_from_rho(real_t rho) const
{
  return np1 * pow( rho/rmd_p, invn ) / h0;
}


/**
\return Specific internal energy \f$ \epsilon = \frac{g-1}{\Gamma} \f$
*/
real_t eos_barotr_gpoly::eps(real_t gm1) const
{
  return h0 * gm1 / gamma + sed0;
}


/**
\return Pressure 
\f$ P = \rho_p \left( h_0 \frac{g-1}{1+n} \right)^{1+n} \f$
*/
real_t eos_barotr_gpoly::press(real_t gm1) const
{
  return rmd_p * pow( h0 * gm1 / np1, np1 );
}

/**
\return Rest mass density 
\f$ \rho = \rho_p \left( h_0 \frac{g-1}{1+n} \right)^n \f$
*/
real_t eos_barotr_gpoly::rho(real_t gm1) const
{
  return rmd_p * pow( h0 * gm1 / np1, n);
}

/**
\return specific enthalpy \f$ h-1 = (g-1) h_0 + \epsilon_0 \f$
*/
real_t eos_barotr_gpoly::hm1(real_t gm1) const
{
  return gm1*h0 + sed0;
}

/**
\return Soundspeed squared \f$ c_s^2 = \frac{h-1}{nh} \f$
*/
real_t eos_barotr_gpoly::csnd(real_t gm1) const
{
  //TODO: check
  real_t h1 = hm1(gm1);
  real_t cs2 = h1/(n*(h1+1));
  return sqrt(cs2);
}


real_t eos_barotr_gpoly::ye(real_t gm1) const
{
  throw std::runtime_error("eos_barotr_gpoly: electron fraction not "
                           "defined for this EOS");
}

real_t eos_barotr_gpoly::rmd_p_from_p_rho_n(real_t p, 
                                      real_t rho, real_t n)
{
  return rho * pow(rho/p, n);
}

real_t eos_barotr_gpoly::eps0_from_p_rho_eps_n(real_t p, 
                                    real_t rho, real_t eps, real_t n)
{
  return eps - n * p / rho;
}


eos_barotr_gpoly 
eos_barotr_gpoly::from_boundary(real_t rho0, real_t eps0,
                                real_t p0, real_t n, real_t rho_max)
{
  return {n, rmd_p_from_p_rho_n(p0, rho0, n), 
          eps0_from_p_rho_eps_n(p0, rho0, eps0, n), rho_max};
}

eos_barotr EOS_Toolkit::make_eos_barotr_gpoly(real_t n, real_t rmd_p, 
                                   real_t sed0, real_t rho_max)
{
  return eos_barotr{std::make_shared<eos_barotr_gpoly>(n, rmd_p,
                                                  sed0, rho_max)};
}
