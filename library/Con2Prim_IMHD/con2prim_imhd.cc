/*! \file con2prim_imhd.cc
\brief Ideal MHD primitive variable recovery algorithm
*/

#include "con2prim_imhd_internals.h"
#include <cassert>
#include <cmath>
#include <limits>
#include "find_roots.h"

using namespace EOS_Toolkit;
using namespace EOS_Toolkit::detail;
using namespace std;


con2prim_mhd::con2prim_mhd(eos_thermal eos_, real_t rho_strict_, 
    bool ye_lenient_, real_t z_lim_, real_t b_lim_, 
    const atmosphere& atmo_, real_t acc_, int max_iter_) 
: eos(std::move(eos_)), rho_strict(rho_strict_), 
  ye_lenient(ye_lenient_), z_lim(z_lim_),
  bsqr_lim(b_lim_*b_lim_), atmo(atmo_), acc(acc_), max_iter(max_iter_)
{
  w_lim = sqrt(1.0 + z_lim*z_lim);
  v_lim = z_lim / w_lim;
}

/**
Requires basic physical constraints
\f[ 
h_0>0, \qquad
r^2 \ge 0, \qquad
(r^l b_l)^2 \ge 0, \qquad
b^2 \ge 0
\f]
The safety margin for the root bracketing needs to satisfy
\f$ 1_m > 1 \f$.
**/
f_upper::f_upper(real_t h0_, real_t rsqr_, 
  real_t rbsqr_, real_t bsqr_)
: h0(h0_), h0sqr(h0_*h0_), rsqr(rsqr_), rbsqr(rbsqr_), bsqr(bsqr_)
{
  assert(h0 > 0);
  assert(rsqr >= 0);
  assert(rbsqr >= 0);
  assert(bsqr >= 0);  
}

/**
Computes 
\f[ x = \frac{y}{y+b^2} = \frac{1}{1 + \mu b^2} \f]
**/
real_t f_upper::x_from_mu(const real_t mu) const
{
  return 1 / (1 + mu * bsqr);
}

/**
Uses the formula
\f[
\bar{r}^2 = x^2 r^2_\perp + r^2_\parallel 
= x \left( r^2 x + \mu \left( x + 1 \right) \left(r^l b_l\right)^2 \right)
\f]
**/
real_t f_upper::rfsqr_from_mu_x(const real_t mu,
  const real_t x) const
{
  return x * (rsqr * x + mu * (x + 1.0) * rbsqr);  
}

/**
Based on the identity 
\f[
hW = \sqrt{h^2 + \bar{r}^2}
\f]
and assuming the minimal value for the enthalpy \f$ h = h_0 \f$
**/
real_t f_upper::new_h0w_from_mu_x(const real_t mu,
  const real_t x) const
{
  return sqrt(h0sqr + rfsqr_from_mu_x(mu,x));  
}


/**
This implements the auxiliary root function as defined in the article
**/
auto f_upper::operator()(const real_t mu) const 
-> std::pair<real_t, real_t>
{
  real_t x     = x_from_mu(mu);
  real_t xsqr  = x*x;
  real_t hw    = new_h0w_from_mu_x(mu, x);
  real_t b     = x * (xsqr * rsqr + mu * (1 + x + xsqr) * rbsqr);
  real_t f     = mu * hw - 1.;
  real_t df    = (h0sqr + b) / hw;
  return {f, df};
}


/**
This computes an initial bracket for the auxiliary root
**/
auto f_upper::initial_bracket() const -> interval<real_t>
{
  real_t mu_min    = 1. / sqrt(h0sqr + rsqr);
  real_t mu0       = 1. / h0;
  real_t rfsqr_min = rfsqr_from_mu_x(mu0, x_from_mu(mu0));
  real_t mu_max    = 1.0 / sqrt(h0sqr + rfsqr_min);
  real_t margin    = 10*std::numeric_limits<real_t>::epsilon();
  mu_max *= 1.0 + margin;
  mu_min *= 1.0 - margin;
  assert(mu_max > mu_min);
  return {mu_min, mu_max};
}


/**
This sets the parameters and EOS defining the root function. We also 
compute the electron fraction limited to the allowed range of the EOS.
Further, we compute an upper limit for the velocity from
\f[
z = \frac{\bar{r}}{h} 
  \le \frac{r}{h} \le \frac{r}{h_0} 
\f]
where we used \f$ \bar{r} \le r \f$ and the minimum enthalpy \f$ h_0 \f$
provided by the EOS. 
**/
froot::froot(const eos_thermal& eos_, real_t valid_ye,
      real_t d_, real_t qtot_, real_t rsqr_, real_t rbsqr_,
      real_t bsqr_, cache& last_ )
: eos(eos_), h0(eos_.minimal_h()), 
  rho_range(eos_.range_rho()), d(d_), qtot(qtot_), rsqr(rsqr_), 
  rbsqr(rbsqr_), bsqr(bsqr_), 
  brosqr(rsqr_ * bsqr_ - rbsqr_), last(last_)
{
  assert(eos.range_ye().contains(valid_ye));
  last.ye    = valid_ye;
  last.calls = 0;
  
  real_t zsqrinf = rsqr / (h0*h0);
  real_t wsqrinf = 1 + zsqrinf;
  winf    = sqrt(wsqrinf);
  vsqrinf = zsqrinf / wsqrinf;
}


/**
Computes 
\f[ x = \frac{y}{y+b^2} = \frac{1}{1 + \mu b^2} \f]
**/
real_t froot::x_from_mu(const real_t mu) const
{
  return 1 / (1 + mu * bsqr);  
}

/**
Uses the formula
\f[
\bar{r}^2 = x^2 r^2_\perp + r^2_\parallel 
= x \left( r^2 x + \mu \left( x + 1 \right) \left(r^l b_l\right)^2 \right)
\f]
**/
real_t froot::rfsqr_from_mu_x(const real_t mu, 
  const real_t x) const
{
  return x * (rsqr * x + mu * (x + 1.0) * rbsqr);  
}

/**
Uses the formula
\f[
\bar{q} = q - \frac{1}{2} \left( b^2 + \mu^2 x^2 b^2 r^2_\perp \right)
\f]
**/
real_t froot::qf_from_mu_x(const real_t mu, 
  const real_t x) const
{
  real_t mux = mu * x;
  return qtot - (bsqr + mux*mux*brosqr) / 2;
}
/**
The formula is written in a way that is accurate also for small velocities.
\f[
\epsilon = W \left( \bar{q} - \mu \bar{r}^2 \right) 
           + v^2 \frac{W^2}{1 + W}
\f]
**/
real_t froot::get_eps_raw(const real_t mu, const real_t qf, 
  const real_t rfsqr, const real_t w) 
{
  return w * (qf - mu * rfsqr*(1.0 - mu * w / (1 + w)));
}

/**
This implements the master root function as defined in the 
article.
**/
real_t froot::operator()(const real_t mu) 
{
  cache& c{last};
  
  c.lmu               = mu;
  c.x                 = x_from_mu(mu);
  const real_t rfsqr  = rfsqr_from_mu_x(mu, c.x);
  const real_t qf     = qf_from_mu_x(mu, c.x);
  c.vsqr              = rfsqr * mu*mu;
  
  
  if (c.vsqr >= vsqrinf) {
    c.vsqr = vsqrinf;
    c.w    = winf;
  } else {
    c.w    = 1 / sqrt(1 - c.vsqr);
  }

  c.rho_raw     = d / c.w;
  c.rho         = rho_range.limit_to(c.rho_raw);

  c.eps_raw     = get_eps_raw(mu, qf, rfsqr, c.w);
  c.eps         = eos.range_eps(c.rho, c.ye).limit_to(c.eps_raw); 

  c.press       = eos.at_rho_eps_ye(c.rho, c.eps, c.ye).press();
  ++c.calls;


  const real_t a        = c.press / (c.rho * (1. + c.eps));

  const real_t h        = (1 + c.eps) * (1 + a);
  
  const real_t hbw_raw  = (1 + a) * (1 + qf - mu * rfsqr);
  const real_t hbw      = max(hbw_raw, h / c.w);     
  const real_t newmu    = 1 / (hbw + rfsqr * mu);
  
  return mu - newmu; 
}


bool froot::stopif(real_t mu, real_t dmu, real_t acc) const
{
  return fabs(dmu) * last.w * last.w < mu * acc;
}


auto froot::initial_bracket(report& errs) const -> interval<real_t>
{
  real_t mu_max {1.0 / h0};

  if (rsqr >= h0*h0) { 
    const int ndigits2{ 36 };
    const real_t margin{ pow(2., 3-ndigits2) };

    f_upper g(h0, rsqr, rbsqr, bsqr);

    ROOTSTAT status;
    mu_max = findroot_using_deriv(g, status, ndigits2, ndigits2+4);

    if (status != ROOTSTAT::SUCCESS) {
      if (status == ROOTSTAT::NOT_CONVERGED) {
        errs.set_prep_root_conv();
      }
      else if (status == ROOTSTAT::NOT_BRACKETED) {
        errs.set_prep_root_bracket();
      }
      return {0, 1.0/h0};
    }
    
    mu_max *= 1. + margin;
    
    assert(g(mu_max).first > 0);
  }
  
  return {0., mu_max};
}


void con2prim_mhd::set_to_nan(prim_vars_mhd& pv, 
                              cons_vars_mhd& cv)
{
  pv.set_to_nan();
  cv.set_to_nan();
}


void con2prim_mhd::operator()(prim_vars_mhd& pv, cons_vars_mhd& cv, 
                               const sm_metric3& g, report& errs) const
{
  errs.iters        = 0;
  errs.adjust_cons  = false;
  errs.set_atmo     = false;
  errs.status       = report::SUCCESS;

  if ((!isfinite(g.vol_elem)) || (g.vol_elem <= 0)) {
    errs.set_invalid_detg(g.vol_elem);
    set_to_nan(pv, cv);
    return;
  } 

  pv.B      = cv.bcons / g.vol_elem;

  const real_t d = cv.dens / g.vol_elem;

  if (d <= atmo.rho_cut) {
    errs.set_atmo_set();
    atmo.set(pv, cv, g);
    return;
  }

  const sm_vec3u bu   = cv.bcons / (g.vol_elem * sqrt(d));
  const sm_vec3l rl   = cv.scon / cv.dens;

  const sm_vec3u ru   = g.raise(rl);
  const real_t rsqr   = ru * rl;
  const real_t rb     = rl * bu;
  const real_t rbsqr  = rb * rb;
  const real_t bsqr   = g.contract(bu, bu);
  const real_t q      = cv.tau / cv.dens;
  const real_t ye0    = cv.tracer_ye / cv.dens;

  if ((!isfinite(d)) || (!isfinite(rsqr))  || (!isfinite(ye0)) ||
      (!isfinite(q)) || (!isfinite(rbsqr)) || (!isfinite(bsqr))) 
  {
    errs.set_nans_in_cons(d, q, rsqr, rbsqr, bsqr, ye0);
    set_to_nan(pv, cv);
    return;
  }     

  if (bsqr > bsqr_lim) 
  {
    errs.set_b_limit(bsqr);
    set_to_nan(pv, cv);
    return; 
  }
  
  const real_t ye = eos.range_ye().limit_to(ye0);


  froot::cache sol{};
  froot f{eos, ye, d, q, rsqr, rbsqr, bsqr, sol}; 

  auto bracket = f.initial_bracket(errs);
  
  if (errs.failed()) {
    set_to_nan(pv, cv);
    return;
  }
  
  rarecase nc(bracket, eos.range_rho(), f);

  if (nc.rho_too_big) 
  {
    errs.set_range_rho(d, d);
    set_to_nan(pv, cv);
    return;
  }

  if (nc.rho_too_small) 
  {
    errs.set_atmo_set();
    atmo.set(pv, cv, g);
    return;
  }

  
  ROOTSTAT status;  
  bracket = findroot_no_deriv(f, nc.bracket, acc, max_iter, status);
  
  errs.iters = sol.calls;
  if (status != ROOTSTAT::SUCCESS) {
    if (status == ROOTSTAT::NOT_CONVERGED) {
      errs.set_root_conv();
    }
    else if (status == ROOTSTAT::NOT_BRACKETED) {
      if (nc.rho_big) { //That's why
        errs.set_range_rho(d, d);
        set_to_nan(pv, cv);
        return;
      }
      if (nc.rho_small) {
        errs.set_atmo_set();
        atmo.set(pv, cv, g);
        return;
      }
      errs.set_root_bracket();
    }
    set_to_nan(pv, cv);
    return;
  }
  assert(bracket.contains(sol.lmu));

  
  if (sol.rho < atmo.rho_cut) {
    errs.set_atmo_set();
    atmo.set(pv, cv, g);
    return;
  }

  auto rgeps = eos.range_eps(sol.rho, sol.ye);
  if (sol.eps_raw > rgeps) {
    errs.adjust_cons = true;
    if (sol.rho >= rho_strict) {
      errs.set_range_eps(sol.eps_raw);
      set_to_nan(pv, cv);
      return;
    }
  }
  else if ( sol.eps_raw < rgeps ) {
    errs.adjust_cons = true;    
  }
    
  
  if (! eos.range_ye().contains(ye0) ) {
    errs.adjust_cons = true;
    if ((!ye_lenient) && (sol.rho >= rho_strict)) {
      errs.set_range_ye(ye0);
      set_to_nan(pv, cv);
      return;
    }
  }

  pv.rho    = sol.rho;
  pv.eps    = sol.eps;
  pv.ye     = sol.ye;
  pv.press  = sol.press;
  pv.vel    = sol.lmu * sol.x  * (ru + (rb * sol.lmu) * bu);
  pv.w_lor  = sol.w;

  real_t sol_v = sqrt(sol.vsqr);
  if (sol_v > v_lim) {
    pv.rho        = d / w_lim;
    if (pv.rho >= rho_strict) {
      errs.set_speed_limit(sol_v);
      set_to_nan(pv, cv);
      return;
    }
    pv.vel       *= v_lim / sol_v;
    pv.w_lor      = w_lim;

    pv.eps = eos.range_eps(pv.rho, pv.ye).limit_to(pv.eps); 

    pv.press       = eos.at_rho_eps_ye(pv.rho, pv.eps, pv.ye).press();
    
    errs.adjust_cons = true;   
  }

  sm_vec3l El = g.cross_product(pv.B, pv.vel);
  pv.E = g.raise(El);

  if (errs.adjust_cons) {
    cv.from_prim(pv, g);
  }
}





f_rare::f_rare(real_t wtarg_, const froot& f_)
: v2targ(1.0 - 1.0/(wtarg_*wtarg_)), f(f_) {}


/**
This implements a root function for finding mu from W_hat
**/
auto f_rare::operator()(const real_t mu) const 
-> std::pair<real_t, real_t>
{
  real_t x     = f.x_from_mu(mu);
  real_t xsqr  = x*x;
  real_t rfsqr = f.rfsqr_from_mu_x(mu,x);
  real_t vsqr  = mu*mu * rfsqr;  
  
  real_t y     = vsqr - v2targ;
  real_t dy = 2.0*mu*x*(xsqr * f.rsqr + mu * (xsqr+x+1.0) * f.rbsqr);
  
  return {y, dy};
}


/**
This checks for the corner case where the density might cross the 
allowed bounds while finding the root of the master function. In
that case, a tighter initial root finding interval is constructed
to guarantee uniqueness.
**/
rarecase::rarecase(const interval<real_t> ibracket,
                   const interval<real_t> rgrho, const froot& f)
{
  
  real_t muc0 = ibracket.min();
  real_t muc1 = ibracket.max();
  const int ndigits = 30;
          
  
  if (f.d > rgrho.max()) {
    real_t wc = f.d / rgrho.max(); 
    if (wc  > f.winf) {
      rho_too_big = true; 
    }
    else {
      f_rare g(wc, f);
      
      if (g(muc1).first <= 0) {
        rho_too_big = true; 
      }
      else {
        if (g(muc0).first < 0) {
          ROOTSTAT status;
          real_t mucc = findroot_using_deriv(g, ibracket, 
                                             status, ndigits, ndigits+2);
          assert(status == ROOTSTAT::SUCCESS); 
          muc0 = max(muc0, mucc);
          rho_big = true;
        }
      }
    }
  }

  if (f.d <  f.winf * rgrho.min()) {
    real_t wc = f.d / rgrho.min(); 
    if (wc < 1) {
      rho_too_small = true; 
    }
    else {
      f_rare g(wc, f);
      if (g(muc0).first >= 0) {
        rho_too_small = true; 
      }
      else {
        if (g(muc1).first > 0) {
          ROOTSTAT status;
          real_t mucc = findroot_using_deriv(g, ibracket, 
                                             status, ndigits, ndigits+2);
          assert(status == ROOTSTAT::SUCCESS); 
          muc1 = min(muc1, mucc);
          rho_small = true;
        }
      }
    }
  }

  bracket = interval<real_t>{muc0, muc1};
}
