/*! \file hydro_prim.h
\brief Class definitions representing primitive hydrodynamic variables.
*/

#ifndef HYDRO_PRIM_H
#define HYDRO_PRIM_H

#include "smtensor.h"

namespace EOS_Toolkit {

/// Structure to represent the primitive hydrodynamic variables.
struct prim_vars {
  real_t rho;      ///< Rest mass density \f$ \rho \f$
  real_t eps;      ///< Specific internal energy \f$ \epsilon \f$
  real_t ye;       ///< Electron fraction \f$ Y_e \f$
  real_t press;    ///< Pressure \f$ P \f$
  sm_vec3u vel;    ///< 3-velocity \f$ v^i \f$
  real_t w_lor;    ///< Lorentz factor \f$ W \f$

  /// Default constructor. Leaves all members uninitialized.
  prim_vars()                            = default;
  
  /// Trivial copy constructor
  prim_vars(const prim_vars&)            = default;
  
  /// Trivial move constructor
  prim_vars(prim_vars&&)                 = default;
  
  /// Trivial copy assignment 
  prim_vars& operator=(const prim_vars&) = default;
  
  /// Trivial move assignment
  prim_vars& operator=(prim_vars&&)      = default;

  /// Construct from single variables.
  prim_vars(real_t rho_, real_t eps_, real_t ye_, 
              real_t press_, sm_vec3u vel_, real_t w_lor_)
  : rho(rho_), eps(eps_), ye(ye_), press(press_), vel(vel_), 
    w_lor(w_lor_) {}
    
  
  /// Convenience method to copy the members into single variables.
  void scatter(real_t& rho_, real_t& eps_, real_t& ye_, real_t& press_, 
               real_t& velx_, real_t& vely_, real_t& velz_, 
               real_t& w_lor_) const;
  
  /// Set all data to NAN
  void set_to_nan();
};


/// Class representing magnetohydrodynamic primitive variables.
struct prim_vars_mhd : public prim_vars {
  sm_vec3u E;       ///< Electric field \f$ E^i \f$
  sm_vec3u B;       ///< Magnetic field \f$ B^i \f$

  /// Default constructor, no initialization.
  prim_vars_mhd() = default;
  
  /// Trivial copy constructor
  prim_vars_mhd(const prim_vars_mhd&)            = default;
  
  /// Trivial move constructor
  prim_vars_mhd(prim_vars_mhd&&)                 = default;
  
  /// Trivial copy assignment 
  prim_vars_mhd& operator=(const prim_vars_mhd&) = default;
  
  /// Trivial move assignment
  prim_vars_mhd& operator=(prim_vars_mhd&&)      = default;

  /// Construct from single variables.
  prim_vars_mhd(real_t rho_, real_t eps_, real_t ye_, 
              real_t press_, sm_vec3u vel_, real_t w_lor_, 
              sm_vec3u E_, sm_vec3u B_)
  : prim_vars(rho_, eps_, ye_, press_, vel_, w_lor_), E(E_), B(B_) {}

  /// Convenience method to copy the members into single variables.
  void scatter(real_t& rho_, real_t& eps_, real_t& ye_, real_t& press_, 
               real_t& velx_, real_t& vely_, real_t& velz_, 
               real_t& w_lor_, real_t& E_x_, real_t& E_y_, real_t& E_z_,
               real_t& B_x_, real_t& B_y_, real_t& B_z_) const;

  /// Set all data to NAN
  void set_to_nan();
};


}
#endif
