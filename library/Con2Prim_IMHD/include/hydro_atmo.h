/*! \file hydro_atmo.h
\brief Class definition representing artificial atmosphere.
*/

#ifndef HYDRO_ATMO_H
#define HYDRO_ATMO_H

#include "hydro_prim.h"
#include "hydro_cons.h"

namespace EOS_Toolkit {


/// Class representing an artificial atmosphere.
struct atmosphere {
  /// Rest mass density of the atmosphere.
  const real_t rho;      
  
  /// Specific internal energy.
  const real_t eps;      

  /// Electron fraction.
  const real_t ye;       

  /// Pressure.
  const real_t press;

  /// Cutoff density below which we set to atmosphere.
  const real_t rho_cut;  
  
  /**\brief Constructor.
  
  @param rho_ Atmosphere mass density
  @param eps_ Atmosphere specific internal energy
  @param ye_  Atmosphere electron fraction
  @param press_ Atmosphere pressure
  @param rho_cut_ Mass density below which to apply atmosphere
  **/
  atmosphere(real_t rho_, real_t eps_, real_t ye_, real_t press_, 
    real_t rho_cut_);
    
  ///Trivial copy constructor  
  atmosphere(const atmosphere&) = default;

  /**\brief Set pure hydro primitive variables to atmosphere.  
  @param pv Primitive variables (hydro)
  **/
  void set(prim_vars& pv) const;

  /**\brief Set pure hydro evolved vars to atmosphere.
  @param cv Evolved variables (hydro)
  @param g  3-metric
  **/
  void set(cons_vars& cv, const sm_metric3& g) const;

  /**\brief Set pure hydro primitive and evolved variables to 
  atmosphere.
  @param pv Primitive variables (hydro)  
  @param cv Evolved variables (hydro)
  @param g  3-metric
  **/   
  void set(prim_vars& pv, cons_vars& cv, 
           const sm_metric3& g) const;

  /** Set ideal MHD primitive variables to atmosphere.
  @param pv Primitive variables (GRMHD)  
  **/
  void set(prim_vars_mhd& pv) const;


  /**\brief Set ideal MHD conserved variables to atmosphere.
  @param cv Evolved variables (GRMHD)
  @param g  3-metric
  **/
  void set(
    cons_vars_mhd& cv,          ///< Conserved vars
    const sm_metric3& g         ///< The 3-metric
  ) const;

  /**\brief Set ideal MHD primitives and conserved variables to 
  atmosphere.
  @param pv Primitive variables (GRMHD)  
  @param cv Evolved variables (GRMHD)
  @param g  3-metric
  **/
  void set(prim_vars_mhd& pv, cons_vars_mhd& cv, 
           const sm_metric3& g) const;

};



}
#endif
