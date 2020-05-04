#ifndef EOS_IDEALGAS_H
#define EOS_IDEALGAS_H

#include "eos_thermal.h"

namespace EOS_Toolkit {

/**\brief Create a classical ideal gas EOS

@return eos_thermal object representing the ideal gas EOS

@param n        Adiabatic index
@param max_eps  Maximum allowed specific internal energy 
                \f$ \epsilon \f$
@param max_rho  Maximum allowed mass density \f$ \rho \f$ 

\rst
The electron fraction is an unused dummy parameter for this EOS,
with allowed range :math:`[0,1]`.
\endrst
**/
eos_thermal make_eos_idealgas(real_t n, real_t max_eps, real_t max_rho);
    
} // namespace EOS_Toolkit

#endif


