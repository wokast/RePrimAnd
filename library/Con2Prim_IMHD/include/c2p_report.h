/*! \file c2p_report.h
\brief Class definitions for status of primitive recovery.
*/

#ifndef CON2PRIM_MHD_ERROR_H
#define CON2PRIM_MHD_ERROR_H

#include <string>
#include "config.h"

namespace EOS_Toolkit {

class con2prim_mhd;
namespace detail {class froot;}

///Struct to represent the outcome of a con2prim call.
/**
This should be checked after con2prim. The various members indicate
success or failure, if the conserved variables were changed, or
artificial atmosphere applied.
*/
struct c2p_mhd_report  {
  public:
  
  ///Failure flags 
  enum err_code {
    SUCCESS,       /*!<Recovery succeeded (possibly after 
                       allowed corrections) */
    INVALID_DETG,  ///<3-metric determinant not positive or NAN/INF
    NANS_IN_CONS,  ///<One or more evolved variables is NAN/INF
    RANGE_RHO,     ///<Mass density outside EOS range
    RANGE_EPS,     ///<Fluid internal energy outside EOS range
    SPEED_LIMIT,   ///<Speed limit exceeded
    RANGE_YE,      ///<Electron fraction outside EOS range
    B_LIMIT,       ///<Magnetization above limit
    ROOT_FAIL_CONV,        /*!<Root finding did not converge (this 
                               should never happen) */ 
    ROOT_FAIL_BRACKET,     /*!<Bracketing of root failed (this 
                               should never happen) */ 
    PREP_ROOT_FAIL_CONV,   /*!<Auxiliary root finding did not 
                               converge (this should never happen) */ 
    PREP_ROOT_FAIL_BRACKET,/*!<Bracketing of auxiliary root 
                               failed (this should never happen) */ 
    ERR_CODE_NOT_SET
    };
  
  ///Default constructor, resulting object invalid.
  c2p_mhd_report() = default;
  
  
  /// Throw an exception with the debug message
  void raise() const;
  
  /**
  @return String with debug information.
  
  This assembles a human readable debug message from the flags
  and numerical values stored internally (for performance reasons, no
  string is ever created before explicitly requested via this function).
  **/
  std::string debug_message() const;
  
  /** 
  @return If the input was invalid according to the error policy.
  **/
  bool failed() const {return status != SUCCESS;}        
  
  
  /// SUCCESS or reason for failure. 
  err_code status{ERR_CODE_NOT_SET};

  /// Whether the conserved variables were adjusted.
  bool adjust_cons{false};     

  /// Whether the artificial atmosphere was enforced.  
  bool set_atmo{false};        

  /// Number of calls to the EOS needed for the root finding.
  unsigned int iters{0};            

  protected:

  /// Set state artificial atmosphere enforced.
  void set_atmo_set();
  
  /// Set error invalid metric determinant. 
  void set_invalid_detg(real_t detg_);
  
  /// Set error NANs in conserved variables. 
  void set_nans_in_cons(real_t dens_, real_t qtot_, real_t rsqr_,
    real_t rbsqr_, real_t bsqr_, real_t ye_);
    
  /// Set error density out of range. 
  void set_range_rho(real_t dens_, real_t rho_);
    
  /// Set error energy out of range. 
  void set_range_eps(real_t eps_);
  
  /// Set error speed limit exceeded. 
  void set_speed_limit(real_t vel_);

  /// Set error limit for b exceeded. 
  void set_b_limit(real_t bsqr_);
  
  /// Set error energy out of range. 
  void set_range_ye(real_t ye_);
  
  /// Set error root finding did not converge. 
  void set_root_conv();

  /// Set error root bracketing faulty. 
  void set_root_bracket();
  
  /// Set error preperatory root finding not converged. 
  void set_prep_root_conv();

  /// Set error preperatory root bracketing faulty. 
  void set_prep_root_bracket();

  private:
  
  /// Conserved density \f$ D \f$.
  real_t dens;

  /// Specific conserved momentum \f$ r^2 = \frac{ S_i S^i}{D^2} \f$.
  real_t rsqr;               

  /// Product \f$ (r^l b_l)^2 \f$.
  real_t rbsqr;              

  /// Specific total energy density \f$ q = \frac{\tau}{D}  \f$.
  real_t qtot;               

  /// Ratio \f$ b^2 = \frac{B^2}{D} \f$.
  real_t bsqr;

  /// Electron fraction \f$ Y_e \f$.
  real_t ye;
  
  /// Rest mass density \f$ \rho \f$.
  real_t rho;
  
  /// Specific internal energy \f$ \epsilon \f$.
  real_t eps;
  
  /// Velocity.
  real_t vel;
  
  ///3-Metric determinant.
  real_t detg; 
  
  
  friend class con2prim_mhd;
  friend class detail::froot;
};

}
#endif
