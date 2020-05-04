/*! \file hydro_cons.cc
\brief Implementation of class for artificial atmosphere.
*/

#include "hydro_atmo.h"
namespace EOS_Toolkit {


atmosphere::atmosphere(real_t rho_, real_t eps_, real_t ye_, real_t press_, real_t rho_cut_)
: rho(rho_), eps(eps_), ye(ye_), press(press_), rho_cut(rho_cut_) {}



void atmosphere::set(prim_vars& pv) const
{
  pv.rho    = rho;
  pv.eps    = eps;
  pv.ye     = ye;
  pv.press  = press;
  pv.vel    = ZERO;
  pv.w_lor  = 1.0;
}

void atmosphere::set(cons_vars& cv, const sm_metric3& g) const
{
  cv.dens       = g.vol_elem * rho; 
  cv.tau        = cv.dens * eps;
  cv.tracer_ye  = cv.dens * ye;
  cv.scon       = ZERO;
}

void atmosphere::set(prim_vars& pv, cons_vars& cv, 
                     const sm_metric3& g) const
{
  set(pv);
  set(cv, g);
}


void atmosphere::set(prim_vars_mhd& pv) const
{
  set((prim_vars&)pv);
  //TODO: decide what to do with EM part
  pv.E = ZERO;
}

void atmosphere::set(cons_vars_mhd& cv, const sm_metric3& g) const
{
  set((cons_vars&) cv, g);
  //TODO: decide what to do with EM part
  cv.tau += 0.5 * g.norm2(cv.bcons) / g.vol_elem; 
}

void atmosphere::set(prim_vars_mhd& pv, cons_vars_mhd& cv, 
                     const sm_metric3& g) const
{
  set(pv);
  set(cv, g);
}


}
