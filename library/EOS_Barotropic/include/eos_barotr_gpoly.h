#ifndef EOS_BAROTR_GPOLY_H
#define EOS_BAROTR_GPOLY_H

#include "eos_barotropic.h"
#include <vector>

namespace EOS_Toolkit {

eos_barotr make_eos_barotr_gpoly(
  real_t n_,                          ///<Adiabatic index \f$ n \f$
  real_t rmd_p_,                      ///<Density scale \f$ \rho_p \f$
  real_t sed0_,                       ///< \f$ \epsilon_0 \f$    
  real_t rho_max_                     ///<Max valid density 
);


}//namespace EOS_Toolkit


#endif

