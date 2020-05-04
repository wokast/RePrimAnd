#ifndef EOS_BAROTR_TABLE_H
#define EOS_BAROTR_TABLE_H

#include "eos_barotropic.h"
#include <vector>

namespace EOS_Toolkit {

/**\brief Create tabulated barotropic EOS

The internally used regularly log-spaced lookup tables are set up
from irregularly spaced sample points provided as vectors, using
cubic monotonic spline interpolation.
 
One can specify whether the EOS is isentropic. 
If the temperature is zero, the EOS must be specified as isentropic 
or an exception is thrown.
Even for isentropic EOS, it is still necessary to provide specific 
energy and pseudo-enthalpy, even though they become redundant. 
Consistency is not checked and in the responsibility of the user.

The validity
range will start at zero up to the largest provided sample point.
Between zero and the lowest sample point, a polytropic EOS will be 
used. 

@param gm1 Samples for pseudo enthalpy \f$ g-1 \f$. Must be strictly 
           increasing.
@param rho Samples for mass density \f$ \rho \f$. Must be strictly 
           increasing.
@param eps Samples for specific energy \f$ \epsilon \f$. Must be 
           strictly increasing and obey \f$ \epsilon > -1 \f$
@param pbr Samples for \f$ P / \rho \ge 0 \f$.
@param cs2 Samples for squared soundspeed. Must obey 
           \f$ 0 \le c_s^2 < 1 \f$
@param temp Samples for temperature or empty vector (equivalent 
            to all-zero). Must obey \f$ T\ge 0 \f$.
@param efrac Electron fraction or empty vector (in which case the EOS
             will not provide electron fraction).
@param isentropic_ Whether the EOS is isentropic, e.g. for degenerate 
                    matter.
@param n_poly_ Polytropic index used to extent EOS to zero density.

@return Generic interface employing tabulated barotropic EOS
**/
eos_barotr make_eos_barotr_table(
  const std::vector<real_t>& gm1,
  const std::vector<real_t>& rho,
  const std::vector<real_t>& eps,
  const std::vector<real_t>& pbr,
  const std::vector<real_t>& cs2,
  const std::vector<real_t>& temp,
  const std::vector<real_t>& efrac, 
  bool isentropic_,            
  real_t n_poly_              
);


}//namespace EOS_Toolkit


#endif

