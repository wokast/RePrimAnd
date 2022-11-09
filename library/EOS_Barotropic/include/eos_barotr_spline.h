#ifndef EOS_BAROTR_SPLINE_H
#define EOS_BAROTR_SPLINE_H

#include "eos_barotropic.h"
#include <vector>
#include <functional>

namespace EOS_Toolkit {

/**\brief Create barotropic EOS based on splines

The internally used regularly log-spaced cardinal cubic splines are 
set up from irregularly spaced sample points provided as vectors, 
using cubic monotonic spline interpolation.
 
One can specify whether the EOS is isentropic. 
If the temperature is zero, the EOS must be specified as isentropic 
or an exception is thrown.
Even for isentropic EOS, it is still necessary to provide specific 
energy and pseudo-enthalpy, even though they become redundant. 
Consistency is not checked and in the responsibility of the user.

The validity
range will start at zero up to the largest provided sample point.
Between zero and a matching point specified by the user, 
a polytropic EOS will be used. Between the matching point and the 
largest sample point, the EOS will be sampled regularly in log-space
with a given number of points per decade. 

The parameters are w.r.t EOS units. The EOS units given in the 
last parameter are stored for bookkeeping.


@param gm1 Samples for pseudo enthalpy \f$ g-1 \f$. Must be strictly 
           increasing.
@param rho Samples for mass density \f$ \rho \f$. Must be strictly 
           increasing.
@param eps Samples for specific energy \f$ \epsilon \f$. Must be 
           strictly increasing and obey \f$ \epsilon > -1 \f$
@param pbr Samples for \f$ P / \rho \ge 0 \f$.
@param csnd Samples for soundspeed. Must obey 
           \f$ 0 \le c_s < 1 \f$
@param temp Samples for temperature or empty vector (equivalent 
            to all-zero). Must obey \f$ T\ge 0 \f$.
@param efrac Electron fraction or empty vector (in which case the EOS
             will not provide electron fraction).
@param isentropic_ Whether the EOS is isentropic, e.g. for degenerate 
                    matter.
@param rg_rho Range of \f$ \rho \f$ where EOS is computed from sample
                 points. Below, a polytrope is used. Range must be 
                 within the provided sample points.
@param n_poly_ Polytropic index used to extent EOS to zero density.
@param units   Unit system (w.r.t. SI) of the EOS
@param pts_per_mag Sample points per magnitude to be used internally
                   for EOS splines
                    
@return Generic interface employing tabulated barotropic EOS
**/
eos_barotr make_eos_barotr_spline(
  const std::vector<real_t>& gm1,
  const std::vector<real_t>& rho,
  const std::vector<real_t>& eps,
  const std::vector<real_t>& press,
  const std::vector<real_t>& csnd,
  const std::vector<real_t>& temp,
  const std::vector<real_t>& efrac, 
  bool isentropic_,            
  interval<real_t> rg_rho,
  real_t n_poly_,
  units units_=units::geom_solar(),
  std::size_t pts_per_mag=200
);


eos_barotr make_eos_barotr_spline(
  std::function<real_t(real_t)> gm1_rho, 
  std::function<real_t(real_t)> rho_gm1, 
  std::function<real_t(real_t)> eps_gm1, 
  std::function<real_t(real_t)> press_gm1, 
  std::function<real_t(real_t)> csnd_gm1, 
  std::function<real_t(real_t)> temp_gm1, 
  std::function<real_t(real_t)> efrac_gm1, 
  bool isentropic, 
  interval<real_t> rg_rho, 
  real_t n_poly,
  units u, 
  std::size_t pts_per_mag=200
);


/**\brief Create barotropic EOS based on splines

This allows to represent existing EOS objects of any type 
as a spline EOS by sampling them
it.

@param eos The EOS to be sampled
@param rg_rho The density range that should be represented by
              interpolating splines
@param n_poly Polytropic index used to extent EOS to zero density.
@param pts_per_mag Sample points per magnitude to be used internally
                   for EOS splines
                    
@return Generic interface employing tabulated barotropic EOS
*/
eos_barotr make_eos_barotr_spline(
  const eos_barotr& eos, 
  interval<real_t> rg_rho, 
  real_t n_poly,
  std::size_t pts_per_mag=200
);

}//namespace EOS_Toolkit


#endif

