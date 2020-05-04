/*! \file hydro_cons.cc
\brief Implementation of classes representing conserved variables.
*/

#include "hydro_cons.h"
#include <limits>

namespace EOS_Toolkit {


/**
\rst
The conserved variables are computed from

.. math::

   D     &= V_c W \rho \\ 
   S_i   &= V_c W^2 \rho h v_i \\
   \tau  &= V_c W^2 \left( v^2 \left( \rho \frac{W}{1 + W} + P \right) 
                           + \rho \epsilon \right) \\
   Y_e^T &= D Y_e 


where :math:`V_c = \sqrt{\det(g_{ij})}` is the 3-metric volume element.
The formula for :math:`\tau` is written in a way which is accurate 
also in the Newtonian limit, i.e. for 
:math:`v^2 \ll 1, \epsilon \ll 1`.

\endrst
**/
void cons_vars::from_prim(const prim_vars& pv, const sm_metric3& g)
{
  const sm_vec3l vl   = g.lower(pv.vel);
  const real_t hrho  = pv.rho * (1.0 + pv.eps) + pv.press;
  const real_t v2    = vl * pv.vel;
	const real_t wl2   = pv.w_lor * pv.w_lor;
  dens                = g.vol_elem * pv.w_lor * pv.rho;
  scon                = (g.vol_elem * wl2 * hrho) * vl;
  tau                 = g.vol_elem * wl2 * (v2*(pv.rho * pv.w_lor / (1.0 + pv.w_lor) 
                          + pv.press) + pv.rho * pv.eps);
  tracer_ye           = dens * pv.ye;
}



void cons_vars::scatter(real_t& dens_, real_t& tau_, 
  real_t& tracer_ye_, real_t& sconx_, real_t& scony_, 
  real_t& sconz_) const
{
  dens_       = dens;
  tau_        = tau;
  tracer_ye_  = tracer_ye;
  sconx_      = scon(0);
  scony_      = scon(1);
  sconz_      = scon(2);
}

void cons_vars::set_to_nan()
{
  dens 
  = tau 
  = tracer_ye 
  = scon(0)  = scon(1)   = scon(2) 
  = std::numeric_limits<real_t>::quiet_NaN();
}
  
void cons_vars_mhd::add_em_part(const sm_vec3u& E, const sm_vec3u& B, 
  const sm_metric3& g)
{
  real_t esqr = g.norm2(E);
  real_t bsqr = g.norm2(B);
  auto ecrossb = g.cross_product(E, B);
  scon += g.vol_elem * ecrossb;
  tau  += g.vol_elem * 0.5 * (esqr + bsqr);
  bcons = g.vol_elem * B;
}


void cons_vars_mhd::from_prim(const prim_vars_mhd& pv, 
  const sm_metric3& g)
{
  cons_vars::from_prim(pv, g);
  add_em_part(pv.E, pv.B, g);
}


void cons_vars_mhd::scatter(real_t& dens_, real_t& tau_, real_t& tracer_ye_, 
               real_t& sconx_, real_t& scony_, real_t& sconz_,
               real_t& bconx_, real_t& bcony_, real_t& bconz_) const
{
  cons_vars::scatter(dens_, tau_, tracer_ye_, sconx_, scony_, sconz_);
  bconx_      = bcons(0);
  bcony_      = bcons(1);
  bconz_      = bcons(2);
}

void cons_vars_mhd::set_to_nan()
{
  cons_vars::set_to_nan();
  bcons(0) = bcons(1) =  bcons(2) 
  = std::numeric_limits<real_t>::quiet_NaN();
}


}

