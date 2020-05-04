#ifndef EOS_BAROTR_POLY_H
#define EOS_BAROTR_POLY_H

#include "eos_barotropic.h"

namespace EOS_Toolkit {

/**\brief Create a polytropic EOS

The valid range always starts at zero density up to user specified 
maximum. If the polytropic index is less than one, the valid
range is reduced avoid superluminal soundspeed, if necessary.

@param n       Polytropic index \f$ n > 0 \f$
@param rmd_p   Polytropic density scale \f$ \rho_p > 0 \f$
@param rho_max Maximum allowed mass density
@return eos_barotr generic interface employing polytropic EOS
**/
eos_barotr make_eos_barotr_poly(real_t n, real_t rmd_p, 
                                real_t rho_max);


}//namespace EOS_Toolkit


#endif

