#include "eos_barotr_pwpoly.h"
#include "eos_barotr_pwpoly_impl.h"
#include <stdexcept>
#include <cmath>
#include <limits>

using namespace std;
using namespace EOS_Toolkit;
using namespace EOS_Toolkit::implementations;

eos_poly_piece::eos_poly_piece(real_t rmd0_, real_t sed0_, 
                 real_t gamma_, real_t rmd_p_)
: rmd0(rmd0_),  gamma(gamma_), rmd_p(rmd_p_) 
{
  n    = 1.0 / (gamma - 1.0);
  np1  = n + 1.0;
  invn = 1.0 / n;
  dsed = sed0_ - n*pow(rmd0/rmd_p, invn);
  gm10 = gm1_from_rho(rmd0);
  p0   = press_from_gm1(gm10);  
}  

/**
\return \f$ g-1 = h-1 = (n+1) \left(\frac{\rho}{\rho_p}\right)^{1/n} 
         + \delta\epsilon\f$
*/
real_t eos_poly_piece::gm1_from_rho(real_t rho) const
{
  return np1 * pow( rho/rmd_p, invn ) + dsed;
}


/**
\return Specific internal energy \f$ \epsilon = 
\frac{g-1-\delta\epsilon}{\Gamma} +\delta\epsilon\f$
*/
real_t eos_poly_piece::eps_from_gm1(
  real_t gm1  ///< \f$ g-1 \f$
) const
{
  return (gm1 - dsed) / gamma + dsed;
}

real_t eos_poly_piece::eps_from_rho(real_t rho) const
{
  return eps_from_gm1(gm1_from_rho(rho));
}



/**
\return Pressure \f$ P = 
\rho_p \left( \frac{g-1-\delta\epsilon}{1+n} \right)^{1+n} \f$
*/
real_t eos_poly_piece::press_from_gm1(
  real_t gm1 ///< \f$ g-1 \f$
) const
{
  return rmd_p * pow((gm1 - dsed) / np1, np1 );
}

/**
\return Rest mass density \f$ \rho 
= \rho_p \left( \frac{g-1-\delta\epsilon}{1+n} \right)^n \f$
*/
real_t eos_poly_piece::rho_from_gm1(
  real_t gm1 ///< \f$ g-1 \f$
) const
{
  return rmd_p * pow( (gm1-dsed) / np1, n);
}

/**
\return specific enthalpy \f$ h-1 = g-1 \f$
*/
real_t eos_poly_piece::hm1_from_gm1(
  real_t gm1 ///< \f$ g-1 \f$
) const
{
  return gm1;
}

/**
\return Soundspeed squared \f$ c_s^2 = \frac{g-1-\delta\epsilon}{ng} \f$
*/
real_t eos_poly_piece::csnd_from_gm1(real_t gm1) const
{
  real_t gm1p  = gm1 - dsed;
  return sqrt(gm1p / ( n * (gm1 + 1.0)));
}

real_t eos_poly_piece::rho_max_save(real_t rho_max) const
{
  if ( n>=1 ) return rho_max;
  const real_t gm1_save{ std::max(gm10, (n + dsed) / (1.0 - n)) };
  const real_t rho_save{ rho_from_gm1(gm1_save) };
  const real_t margin { 10*std::numeric_limits<real_t>::epsilon() };
  
  return std::min(rho_max, rho_save * (1.0 - margin));
}

bool eos_poly_piece::rho_save_up_to(real_t rho) const
{
  if (n >= 1) return true;
  return csnd_from_gm1(gm1_from_rho(rho)) < 1;
}


eos_barotr_pwpoly::eos_barotr_pwpoly(real_t rmdp0, 
  const vector<real_t>& segm_bound, 
  const  vector<real_t>& segm_gamma,
  real_t rho_max_)
{
  if (rho_max_ <= 0) {
    throw runtime_error("eos_barotr_pwpoly: maximum density must be "
                        "strictly positive");
  }
  if (segm_bound.size() != segm_gamma.size()) {
    throw runtime_error("eos_barotr_pwpoly: vector sizes mismatch.");
  }
  if (segm_bound.empty()) {
    throw runtime_error("eos_barotr_pwpoly: need at least one "
                        "segment.");
  }
  if (segm_bound[0] != 0) {
    throw runtime_error("eos_barotr_pwpoly: First segment has to "
                        "start at zero density.");
  }
  
  for (std::size_t i = 1; i < segm_bound.size(); ++i) 
  {
    if (segm_bound[i] <= segm_bound[i-1]) {
      throw runtime_error("eos_barotr_pwpoly: segment boundary "
                          "densities not strictly increasing.");      
    }
  }
  
  auto ibnd = segm_bound.begin();
  auto iga  = segm_gamma.begin();

  segments.emplace_back(0, 0, *iga, rmdp0);

  while ((++ibnd != segm_bound.end()) 
         && (++iga != segm_gamma.end())) 
  {
    const eos_poly_piece& lseg = segments.back();
    
    if (*ibnd > rho_max_) break;
    
    if (!lseg.rho_save_up_to(*ibnd)) {
      rho_max_ = lseg.rho_max_save(rho_max_);
      break; 
    }

    real_t sedc = lseg.eps_from_rho(*ibnd);
    real_t np   = 1.0 / (*iga - 1.0);
    real_t et   = np / lseg.n;
    real_t rmdp = pow(lseg.rmd_p, et) * pow(*ibnd, 1.0-et);

    segments.emplace_back(*ibnd, sedc, *iga, rmdp);
  }
  
  
  rgrho = {0, rho_max_};
  rggm1 = {0, gm1_from_rho(rho_max_)};
}


const eos_poly_piece& 
eos_barotr_pwpoly::segment_for_rho(real_t rho) const
{
  auto i = segments.rbegin();
  while (i->rmd0 > rho) {
    if (++i == segments.rend()) return segments[0];
  }
  return *i;
}


const eos_poly_piece& 
eos_barotr_pwpoly::segment_for_gm1(real_t gm1) const
{
  auto i = segments.rbegin();
  while (i->gm10 > gm1) {
    if (++i == segments.rend()) return segments[0];
  }
  return *i;
}
  

/**
\return \f$ g-1 \f$ from polytropic piece at given density.
*/
real_t eos_barotr_pwpoly::gm1_from_rho(real_t rho) const
{
  return segment_for_rho(rho).gm1_from_rho(rho);
}



/**
\return Specific internal energy \f$ \epsilon \f$
from polytropic piece at given \f$ g-1 \f$.
*/
real_t eos_barotr_pwpoly::eps(
  real_t gm1  ///< \f$ g-1 \f$
) const
{
  return segment_for_gm1(gm1).eps_from_gm1(gm1);
}


/**
\return Pressure \f$ P \f$
from polytropic piece at given \f$ g-1 \f$.
*/
real_t eos_barotr_pwpoly::press(
  real_t gm1 ///< \f$ g-1 \f$
) const
{
  return segment_for_gm1(gm1).press_from_gm1(gm1);
}

/**
\return Rest mass density \f$ \rho \f$
from polytropic piece at given \f$ g-1 \f$.
*/
real_t eos_barotr_pwpoly::rho(
  real_t gm1 ///< \f$ g-1 \f$
) const
{
  return segment_for_gm1(gm1).rho_from_gm1(gm1);
}

/**
\return specific enthalpy \f$ h-1 = g-1 \f$.
*/
real_t eos_barotr_pwpoly::hm1(
  real_t gm1 ///< \f$ g-1 \f$
) const
{
  return gm1;
}

/**
\return Soundspeed \f$ c_s  \f$
from polytropic piece at given \f$ g-1 \f$.
*/
real_t eos_barotr_pwpoly::csnd(real_t gm1) const
{
  return segment_for_gm1(gm1).csnd_from_gm1(gm1);
}

real_t eos_barotr_pwpoly::ye(real_t gm1) const
{
  throw std::runtime_error("eos_barotr_pwpoly: electron fraction not "
                           "defined for this EOS");
}
  
  
EOS_Toolkit::eos_barotr 
EOS_Toolkit::make_eos_barotr_pwpoly(real_t rmdp0,  
  const std::vector<real_t>& segm_bound,  
  const std::vector<real_t>& segm_gamma,  real_t rho_max)
{
  return eos_barotr{std::make_shared<eos_barotr_pwpoly>(rmdp0, 
                                segm_bound, segm_gamma, rho_max)};
}

