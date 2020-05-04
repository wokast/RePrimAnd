/*! \file hydro_cons.h
\brief Class definitions representing conserved hydrodynamic variables.
*/

#ifndef HYDRO_CONS_H
#define HYDRO_CONS_H

#include "hydro_prim.h"

namespace EOS_Toolkit {

/// Structure to represent conserved variables for pure hydrodynamics.
struct cons_vars {
  /// Conserved density \f$ D \f$
  real_t dens;       
  
  /// Conserved energy \f$ \tau \f$.
  real_t tau;        
  
  /// Conserved tracer \f$ Y_e^T \f$ for the electron fraction.
  real_t tracer_ye;  

  /// Conserved momentum \f$ S_i \f$
  sm_vec3l scon;      

  /// Default constructor, no initialization.
  cons_vars()                 = default;
  
  /// Trivial copy constructor
  cons_vars(const cons_vars&) = default;
  
  /// Trivial move constructor
  cons_vars(cons_vars&&)      = default;
  
  /// Trivial copy assignment 
  cons_vars& operator=(cons_vars&&)      = default;
  
  /// Trivial move assignment
  cons_vars& operator=(const cons_vars&) = default;
  
  /// Construct from single variables.
  cons_vars(real_t dens_, real_t tau_, real_t tracer_ye_, 
                      sm_vec3l scon_)
  : dens{dens_}, tau{tau_}, tracer_ye{tracer_ye_}, scon{scon_} {}

  /// Compute conserved variables from primitives and 3-metric
  void from_prim(const prim_vars& pv, const sm_metric3& g);

  ///Convenience method to copy all members into single variables
  void scatter(real_t& dens_, real_t& tau_, real_t& tracer_ye_, 
               real_t& sconx_, real_t& scony_, real_t& sconz_) const;
 
  ///Set all data to NAN
  void set_to_nan();
};


/// Class to represent ideal MHD conserved variables. 
struct cons_vars_mhd : public cons_vars {
  /// Densitized magnetic field
  sm_vec3u bcons;   

  /// Default constructor, no initialization.
  cons_vars_mhd()                                = default;
  
  /// Trivial copy constructor
  cons_vars_mhd(const cons_vars_mhd&)            = default;
  
  /// Trivial move constructor
  cons_vars_mhd(cons_vars_mhd&&)                 = default;
  
  /// Trivial copy assignment 
  cons_vars_mhd& operator=(cons_vars_mhd&&)      = default;
  
  /// Trivial move assignment
  cons_vars_mhd& operator=(const cons_vars_mhd&) = default;


  /// Construct from single variables.
  cons_vars_mhd(real_t dens_, real_t tau_, real_t tracer_ye_, 
                          sm_vec3l scon_, sm_vec3u bcons_)
  : cons_vars{dens_, tau_, tracer_ye_, scon_}, bcons{bcons_} {}


  /// Compute conserved variables from primitives and 3-metric
  void from_prim(const prim_vars_mhd& pv, const sm_metric3& g);


  ///Convenience method to copy all members into single variables
  void scatter(real_t& dens_, real_t& tau_, real_t& tracer_ye_, 
               real_t& sconx_, real_t& scony_, real_t& sconz_,
               real_t& bconx_, real_t& bcony_, real_t& bconz_) const;

                   
  ///Set all data to NAN
  void set_to_nan();

  private:
  
  /// Adds the electromagnetic part.
  void add_em_part(const sm_vec3u& E, const sm_vec3u& B, 
                   const sm_metric3& g);
};



}
#endif
