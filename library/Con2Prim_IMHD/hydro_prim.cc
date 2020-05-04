/** \file hydro_prim.cc
\brief Implementation of classes for collecting primitive variables
**/

#include "hydro_prim.h"
#include <limits>


namespace EOS_Toolkit {



void prim_vars::scatter(real_t& rho_, real_t& eps_, real_t& ye_, 
       real_t& press_, real_t& velx_, real_t& vely_, real_t& velz_, 
       real_t& w_lor_) const
{
  rho_    = rho; 
  eps_    = eps; 
  ye_     = ye; 
  press_  = press;
  velx_   = vel(0); 
  vely_   = vel(1); 
  velz_   = vel(2);
  w_lor_  = w_lor;
}

void prim_vars::set_to_nan()
{
  rho 
  = eps 
  = ye 
  = press 
  = vel(0) = vel(1) = vel(2) 
  = w_lor 
  = std::numeric_limits<real_t>::quiet_NaN();
}


void prim_vars_mhd::scatter(real_t& rho_, real_t& eps_, real_t& ye_, 
      real_t& press_, real_t& velx_, real_t& vely_, real_t& velz_, 
      real_t& w_lor_, real_t& E_x_, real_t& E_y_, real_t& E_z_,
      real_t& B_x_, real_t& B_y_, real_t& B_z_) const
{
  prim_vars::scatter(rho_, eps_, ye_, press_, 
                     velx_, vely_, velz_, w_lor_);
  E_x_ = E(0);
  E_y_ = E(1);
  E_z_ = E(2);    
  B_x_ = B(0);
  B_y_ = B(1);
  B_z_ = B(2);    
}

void prim_vars_mhd::set_to_nan()
{
  prim_vars::set_to_nan();
  E(0) = E(1) = E(2)
  = B(0) = B(1) = B(2) 
  = std::numeric_limits<real_t>::quiet_NaN();
}

} // namespace EOS_Toolkit
