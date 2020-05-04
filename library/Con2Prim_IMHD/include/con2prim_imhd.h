/*! \file con2prim_imhd.h
\brief Class definitions for primitive recovery interface.
*/

#ifndef CON2PRIM_IMHD_H
#define CON2PRIM_IMHD_H

#include "smtensor.h"
#include "hydro_prim.h"
#include "hydro_cons.h"
#include "hydro_atmo.h"
#include "c2p_report.h"
#include "eos_thermal.h"
#include <string>

namespace EOS_Toolkit {



/**\brief Class representing conservative to primitive conversion 
          for ideal MHD

This function object stores EOS, required accuracy, error policies,
and artificial atmosphere settings.
**/
class con2prim_mhd {
  public:
  
  using range = eos_thermal::range;
  using report = c2p_mhd_report;
  

  /**\brief Constructor
  @param eos_          The EOS
  @param  rho_strict_  Density above which most corrections 
                       are forbidden (strict regime)
  @param ye_lenient_   Whether to allow restricting the electron 
                       fraction to valid range also in the strict 
                       regime
  @param z_lim_        Speed limit in terms of \f$ z = W v \f$
  @param b_lim_        Fail for magnetization \f$ b>b_\mathrm{lim} \f$
  @param atmo_         Specifies artificial atmosphere
  @param acc_          Required accuracy \f$ \Delta \f$ (see article).
  @param max_iter_     Maximum allowed iterations for root finding.
  **/
  con2prim_mhd(eos_thermal eos_, real_t rho_strict_, bool ye_lenient_,         
    real_t z_lim_, real_t b_lim_, const atmosphere& atmo_,  
    real_t acc_, int max_iter_);

  /**\brief Convert from conserved to primitive variables
  
  @param pv  Recovered primitive variables will be stored here
  @param cv  Evolved variables for which to recover primitives. If
             corrections are applied, this contains the corrected values
             after the call.
  @param g   The 3-metric
  @param errs Reports the outcome (validity, corrections, etc).
  
  \rst
  After the call, `pv` and `cv` are fully consistent. 
  If the input was 
  invalid and the error policy did not allow to correct it, 
  `pv` and `cv` members are all set to NAN, and 
  `rep.failed() == true`. If the evolved variables were corrected, 
  `rep.adjust_cons == true`. If artificial atmosphere was applied, 
  `rep.set_atmo == true` in addition.  
  \endrst
  **/
  void operator()(prim_vars_mhd& pv, cons_vars_mhd& cv, 
                  const sm_metric3& g, report& errs) const;

  /// Get prescribed accuracy
  real_t get_acc() const {return acc;}

  /// Get prescribed limit on z
  real_t get_z_lim() const {return z_lim;}

  /// Get prescribed limit on v
  real_t get_v_lim() const {return v_lim;}

  /// Get prescribed limit on b
  real_t get_b_lim() const {return std::sqrt(bsqr_lim);}
  
  /// Get prescribed atmosphere 
  const atmosphere& get_atmo() const {return atmo;}
  
  private:
  
  const eos_thermal eos;         
  const real_t rho_strict;       
  const bool ye_lenient;         
  real_t v_lim;                    
  real_t w_lim;                    
  const real_t z_lim;            
  const real_t bsqr_lim;         
  const atmosphere atmo;         
  const real_t acc;              
  const int max_iter;            


  /// Set primitives and conserved to NaN
  static void set_to_nan(prim_vars_mhd& pv, cons_vars_mhd& cv);

};

}


#endif

