#ifndef EOS_IDEALGAS_H
#define EOS_IDEALGAS_H

#include "eos_thermal.h"

namespace EOS_Toolkit {

/**\brief Create a classical ideal gas EOS

The parameters are w.r.t EOS units. The EOS units given in the 
last parameter are stored for bookkeeping.

@return eos_thermal object representing the ideal gas EOS

@param n        Adiabatic index
@param max_eps  Maximum allowed specific internal energy 
                \f$ \epsilon \f$
@param max_rho  Maximum allowed mass density \f$ \rho \f$ 
@param eos_units Unit system of the EOS. Assumed geometric.

\rst
The electron fraction is an unused dummy parameter for this EOS,
with allowed range :math:`[0,1]`.
\endrst
**/
eos_thermal make_eos_idealgas(real_t n, real_t max_eps, real_t max_rho,
                              units eos_units=units::geom_solar());
    
} // namespace EOS_Toolkit

#endif


