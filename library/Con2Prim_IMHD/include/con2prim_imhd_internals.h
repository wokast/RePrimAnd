/*! \file con2prim_imhd_internals.h
\brief Internal stuff for primitive recovery algorithm.
*/

#ifndef CON2PRIM_IMHD_IMPL_H
#define CON2PRIM_IMHD_IMPL_H

#include "con2prim_imhd.h"

namespace EOS_Toolkit {
namespace detail {

class rarecase;
class f_rare;

/// Function object representing the root function.
/** This contains all the fixed parameters defining the function.
    It also remembers intermediate results from the last evaluation,
    which we do not want to recompute for performance reasons.
    In particular, it remembers if specific internal energy and density
    were in the validity region of the EOS, or if they needed to be limited.
**/
class froot {
  using range   = eos_thermal::range;
  using report  = c2p_mhd_report;

  const eos_thermal eos;    ///< The EOS.
  const real_t h0;         ///< Lower bound for enthalpy, \f$ h_0 \f$
  const range rho_range;    ///< Valid density interval of the EOS. 
  const real_t d;          ///< \f$ d = \frac{D}{\sqrt{\det(g_{ij})}} \f$
  const real_t qtot;       ///< \f$ q = \frac{\tau}{D}  \f$
  const real_t rsqr;       ///< \f$ r^2 = \frac{ S_i S^i}{D^2} \f$
  const real_t rbsqr;      ///< \f$ (r^l b_l)^2 \f$
  const real_t bsqr;       ///< \f$ b^2 = \frac{B^2}{D} \f$
  const real_t brosqr;     ///< \f$ b^2 r^2_\perp = b^2 r^2 - (r^lb_l)^2 \f$
  real_t winf;             ///< Upper bound for Lorentz factor
  real_t vsqrinf;          ///< Upper bound for squared velocity 
  
  real_t x_from_mu(real_t mu) const;
  
  ///Computes fluid momentum from total one, for a given \f$ \mu, x \f$.
  real_t rfsqr_from_mu_x(real_t mu, real_t x) const;
  
  ///Computes fluid energy from total one, for a given \f$ \mu, x \f$.
  real_t qf_from_mu_x(real_t mu, real_t x) const;
  
  ///Computes specific energy from conserved fluid variables and velocity. 
  static real_t get_eps_raw(real_t mu,  real_t qf, 
                            real_t rfsqr, real_t w);  
  
  friend class rarecase;
  friend class f_rare;
     
  public:
  
  using value_t = real_t; 


  
  ///Store intermediate results obtained when evaluating function.
  struct cache {
    real_t ye;       ///< Electron fraction in valid range.
    real_t lmu;      ///< \f$ \mu = \frac{1}{h W} \f$.  
    real_t x;         ///< \f$ \frac{1}{1 + \mu b^2} \f$.
    real_t rho;      ///< Rest mass density \f$ \rho \f$.
    real_t rho_raw;  ///< Rest mass density not limited to EOS range.
    real_t eps;      ///< Specific internal energy, limited to EOS range.
    real_t eps_raw;  ///< Specific internal energy, not limited.  
    real_t press;    ///< Pressure \f$ P \f$.
    real_t vsqr;     ///< Squared 3-velocity \f$ v^2 \f$.  
    real_t w;         ///< Lorentz factor \f$ W \f$.
    unsigned int calls;
  };
  
  /// Constructor
  froot(
    const eos_thermal& eos_,       ///< The EOS
    real_t valid_ye,               ///< Electron fraction
    real_t d_,                     ///< \f$ d = \frac{D}{\sqrt{\det(g_{ij})}} \f$
    real_t qtot_,                  ///< \f$ q = \frac{\tau}{D} + EM\f$
    real_t rsqr_,                  ///< \f$ r^2 = \frac{ S_i S^i}{D^2} \f$
    real_t rbsqr_,                 ///< \f$ (r^l b_l)^2 \f$
    real_t bsqr_,                  ///< \f$ b^2 = \frac{B^2}{D} \f$
    cache& last_                   ///< cache for intermediate results
  );

  froot(const froot&)            = default;
  froot(froot&&)                 = default;


  /// The root function
  real_t operator()(real_t mu);

  /// The convergence criterion for root finding.
  bool stopif(real_t mu, real_t dmu, real_t acc) const;

  /// Initial guess for root finding
  auto initial_bracket(report& errs) const -> interval<real_t>;
  
  private:
  
  cache& last;


};

///Class representing the auxiliary root function 
class f_upper {
  public:
  using value_t = real_t;
  
  /// Constructor
  f_upper(
    real_t h0_,         ///< Lower bound for enthalpy
    real_t rsqr_,       ///< \f$ r^2 = \frac{ S_i S^i}{D^2} \f$
    real_t rbsqr_,      ///< \f$ (r^l b_l)^2 \f$
    real_t bsqr_        ///< \f$ b^2 = \frac{B^2}{D} \f$
  );

  
  /// The function and first derivative
  auto operator()(real_t mu) const -> std::pair<real_t,real_t>;


  /// Initial bracket for root finding
  auto initial_bracket() const -> interval<real_t>;
  
  private:
  
  auto x_from_mu(real_t mu) const -> real_t;
  
  auto rfsqr_from_mu_x(real_t mu, real_t x) const -> real_t;
  
  /// Compute new \f$ h_0 W \f$ from initial guess.
  auto new_h0w_from_mu_x(real_t mu, real_t x) const -> real_t;
  
  const real_t h0;      ///< Lower bound for enthalpy, \f$ h_0 \f$
  const real_t h0sqr;   ///< \f$ h_0^2 \f$
  const real_t rsqr;    ///< Fixed parameter \f$ r^2 = \frac{ S_i S^i}{D^2} \f$
  const real_t rbsqr;   ///< Fixed parameter \f$ (r^l b_l)^2 \f$
  const real_t bsqr;    ///< Fixed parameter \f$ b^2 = \frac{B^2}{D} \f$
};




///Root function used by rarecase class   
class f_rare {
  const real_t v2targ;  ///< Target squared velocity   
  const froot& f;

  public:
  using value_t = real_t;
  
  f_rare(
    real_t wtarg_,      ///< Target Lorentz factor
    const froot& f_
  );

  auto operator()(real_t mu) const -> std::pair<real_t,real_t>;  
};

///Class for handling rare corner case 
class rarecase {  
  public:
  rarecase(
    const interval<real_t> bracket,   ///< Initial master root bracket 
    const interval<real_t> rgrho,     ///< Allowed density range
    const froot& f                    ///< Master root function
  );

  /// Root bracket on which solution is unique
  interval<real_t> bracket;  
  
  bool rho_too_big { false };    ///< Density definitely too large
  bool rho_big { false };        ///< Possibly too large
  bool rho_too_small { false };  ///< Density definitely too small
  bool rho_small { false };      ///< Possibly too small
};

}
}

#endif
