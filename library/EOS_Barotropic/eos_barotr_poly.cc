#include "eos_barotr_poly.h"
#include "eos_barotr_poly_impl.h"
#include <cmath>
#include <limits>
#include <stdexcept>

using namespace std;
using namespace EOS_Toolkit;
using namespace EOS_Toolkit::implementations;

void eos_barotr_poly::init(real_t n_, real_t rmd_p_, real_t rho_max_) 
{
  if (n_<0) {
    throw std::range_error("eos_barotr_poly: "
                           "negative polytropic index"); 
  }
  if (rho_max_ <= 0) {
    throw runtime_error("eos_barotr_poly: maximum density must be "
                        "strictly positive");
  }
  n     = n_;
  rmd_p = rmd_p_;
  np1   = n + 1; 
  gamma = 1.0 + 1.0 / n; 
  invn  = 1.0 / n; 
  
  real_t gm1_max = gm1_from_rho(rho_max_);
  if (n < 1) {
    const real_t gm1_c { n / (1.0 - n) };
    const real_t margin { 10*std::numeric_limits<real_t>::epsilon() };
    gm1_max  = std::min(gm1_max, gm1_c * (1.0 - margin));
    rho_max_ = rho(gm1_max);
  }
  
  rgrho = {0, rho_max_};
  rggm1 = {0, gm1_max};
}


eos_barotr_poly::eos_barotr_poly(real_t n_, real_t rmd_p_, 
                             real_t rho_max_) 
{
  init(n_, rmd_p_, rho_max_);
}

/**
Construct polytrope by specifying \f$ \epsilon \f$ and \f$ P \f$
at a given density \f$ \rho \f$, using the formulas
\f{eqnarray*}{ 
  n      &=& \frac{\epsilon \rho}{P}  \\
  \rho_p &=& \rho \left(\frac{\rho}{P} \right)^n
\f}
*/
eos_barotr_poly::eos_barotr_poly(real_t rmd_m, real_t sed_m, real_t p_m, 
                             real_t rho_max_) 
{
  real_t n_     = sed_m * rmd_m / p_m;
  real_t rmd_p_ = rmd_m * pow(rmd_m / p_m, n_);
  init(n_, rmd_p_, rho_max_);
}

/**
\return \f$ g-1 = h-1 = (n+1) \left(\frac{\rho}{\rho_p}\right)^{1/n} \f$
*/
real_t eos_barotr_poly::gm1_from_rho(real_t rho) const
{
  return np1 * pow( rho/rmd_p, invn );
}


/**
\return Specific internal energy \f$ \epsilon = \frac{g-1}{\Gamma} \f$
*/
real_t eos_barotr_poly::eps(real_t gm1) const
{
  return gm1 / gamma;
}


/**
\return Pressure \f$ P = \rho_p \left( \frac{g-1}{1+n} \right)^{1+n} \f$
*/
real_t eos_barotr_poly::press(real_t gm1) const
{
  return rmd_p * pow( gm1 / np1, np1 );
}

/**
\return Rest mass density \f$ \rho = \rho_p \left( \frac{g-1}{1+n} \right)^n \f$
*/
real_t eos_barotr_poly::rho(real_t gm1) const
{
  return rmd_p * pow( gm1 / np1, n);
}

/**
\return specific enthalpy \f$ h-1 = g-1 \f$
*/
real_t eos_barotr_poly::hm1(real_t gm1) const
{
  return gm1;
}

/**
\return Soundspeed squared \f$ c_s^2 = \frac{g-1}{ng} \f$
*/
real_t eos_barotr_poly::csnd(real_t gm1) const
{
  real_t cs2 = gm1/(n*(gm1+1));
  return sqrt(cs2);
}


real_t eos_barotr_poly::ye(real_t gm1) const
{
  throw std::runtime_error("eos_barotr_poly: electron fraction not "
                           "defined for this EOS");
}


  
eos_barotr EOS_Toolkit::make_eos_barotr_poly(real_t n, real_t rmd_p, 
                                             real_t rho_max)
{
  return eos_barotr{std::make_shared<eos_barotr_poly>(n, 
                                           rmd_p, rho_max)};
}
