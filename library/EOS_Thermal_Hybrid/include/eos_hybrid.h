#ifndef EOSHYBRID_H
#define EOSHYBRID_H
#include "eos_barotropic.h"
#include "eos_thermal.h"

namespace EOS_Toolkit {

/**\brief Create Hybrid EOS
@return eos_thermal object representing the hybrid EOS

@param eos_c    The cold (barotropic) EOS
@param gamma_th The adiabatic exponent \f$ \Gamma_\mathrm{th} \f$
                for the thermal part
@param eps_max  Maximum allowed specific internal energy
                \f$ \epsilon \f$
@param rho_max  Maximum allowed mass density \f$ \rho \f$

\rst
The electron fraction is an unused dummy parameter for this EOS,
with allowed range :math:`[0,1]`.
\endrst
**/
eos_thermal make_eos_hybrid(eos_barotr eos_c,  
             real_t gamma_th, real_t eps_max, real_t rho_max);


} // namespace EOS_Toolkit

#endif


