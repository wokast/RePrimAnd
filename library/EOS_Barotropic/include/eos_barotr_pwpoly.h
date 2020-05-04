#ifndef EOS_BAROTR_PWPOLY_H
#define EOS_BAROTR_PWPOLY_H

#include "eos_barotropic.h"
#include <vector>

namespace EOS_Toolkit {

/**\brief Create piecewise polytropic EOS

The valid range always starts at zero density up to user specified 
maximum. If the polytropic index of some segments is less than one, 
the valid range is reduced, if necessary, to avoid superluminal 
soundspeed.

@param rmdp0 Polytropic density scale \f$\rho_{p,0} > 0\f$ 
             of first segment
@param segm_bound Segment boundaries mass densities
@param segm_gamma Polytropic exponents \f$ \Gamma_i > 1 \f$ for all
                  segments
@param rho_max_   Maximum allowed mass density 
@return eos_barotr generic interface employing piecewise polytropic EOS
**/
eos_barotr make_eos_barotr_pwpoly(real_t rmdp0,
  const std::vector<real_t>& segm_bound,
  const std::vector<real_t>& segm_gamma,
  real_t rho_max_   
);


}//namespace EOS_Toolkit


#endif

