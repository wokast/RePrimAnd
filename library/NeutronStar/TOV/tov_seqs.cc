#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>
#include <stdexcept>
#include <cmath>
#include <memory>
#include <algorithm>
#include <limits>
#include "interpol.h"
#include "intervals.h"
#include "spherical_stars.h"
#include "tov_seqs_impl.h"


using namespace std;


namespace EOS_Toolkit {

auto find_rhoc_tov_max_mass(eos_barotr eos, 
                       const real_t rhobr0, const real_t rhobr1,
                       const int bits, const real_t acc, 
                       unsigned int max_steps)
-> real_t
{
  const real_t rho0 { eos.range_rho().limit_to(rhobr0) };
  const real_t rho1 { eos.range_rho().limit_to(rhobr1) };
  
  tov_acc_simple accs{acc};
  
  auto fmin = [&] (real_t rho_c) {
    const auto tov{ 
      get_tov_star_properties(eos, rho_c, accs, false, false) 
    };
    return -tov.grav_mass();
  };
  
  boost::uintmax_t iters{ max_steps };
  const auto r{ 
    boost::math::tools::brent_find_minima(fmin, rho0, rho1, 
                                          bits, iters)
  };
  
  if (iters >= max_steps)
  {
    throw std::runtime_error("TOV maximum mass not found");
  }
  
  return r.first;
  
}



auto find_rhoc_tov_of_mass(eos_barotr eos, real_t mg, 
                        const real_t rhobr0, const real_t rhobr1,
                        real_t acc, unsigned int max_steps)
-> real_t


{
  tov_acc_simple accs{acc};
  
  auto froot = [&] (real_t rho_c) {
    auto tov = get_tov_star_properties(eos, rho_c, accs, false, false);
    return tov.grav_mass() - mg;
  };

  auto stopif = [&] (real_t l, real_t r) {
    return fabs(l-r) < fabs(l+r) * acc;
  };
  
  boost::uintmax_t iters{ max_steps };
    
  auto res = boost::math::tools::toms748_solve(froot, rhobr0, rhobr1, 
                                               stopif, iters);
  
  if (iters >= max_steps)
  {
    throw std::runtime_error("TOV model with mass: root finding failed");
  }
  
  return res.first;
}

namespace detail {

star_seq_impl::star_seq_impl(spline_t mg_gm1_, spline_t mb_gm1_, 
                spline_t rc_gm1_, spline_t mi_gm1_, 
                spline_t lt_gm1_, range_t rg_gm1_,
                units u_)
: mg_gm1{std::move(mg_gm1_)}, mb_gm1{std::move(mb_gm1_)}, 
  rc_gm1{std::move(rc_gm1_)}, mi_gm1{std::move(mi_gm1_)},
  lt_gm1{std::move(lt_gm1_)}, rg_gm1{std::move(rg_gm1_)}, 
  u{u_}
{
  
  if (mg_gm1.range_x().min() <= 0)
  {
    throw std::runtime_error("Attempt to create star sequence "
                              "with invalid pseudo enthalpy");
  }  
  if (mg_gm1.range_y().min() <= 0)
  {
    throw std::runtime_error("Attempt to create star sequence "
                              "with negative grav. mass");
  }  
  if (mb_gm1.range_y().min() <= 0)
  {
    throw std::runtime_error("Attempt to create star sequence "
                              "with negative baryonic mass");
  }  
  if (rc_gm1.range_y().min() <= 0)
  {
    throw std::runtime_error("Attempt to create star sequence "
                      "with negative proper circumferential radius");
  }  
}

auto star_seq_impl::from_vector(std::vector<real_t> mg, 
       std::vector<real_t> mb, std::vector<real_t> rc, 
       std::vector<real_t> mi, std::vector<real_t> lt, 
       range_t rg_gm1, units u)
-> std::shared_ptr<star_seq_impl>
{

  auto mg_gm1 = make_interpol_regspl(std::move(mg), rg_gm1);
  auto mb_gm1 = make_interpol_regspl(std::move(mb), rg_gm1);
  auto rc_gm1 = make_interpol_regspl(std::move(rc), rg_gm1);
  auto mi_gm1 = make_interpol_regspl(std::move(mi), rg_gm1);
  auto lt_gm1 = make_interpol_regspl(std::move(lt), rg_gm1);
  
  
  return make_shared<star_seq_impl>(mg_gm1, mb_gm1, rc_gm1,
                                    mi_gm1, lt_gm1, rg_gm1, u);
}


auto star_seq_impl::grav_mass_from_center_gm1(real_t gm1c) 
const -> real_t
{
  return mg_gm1(gm1c);
}

auto star_seq_impl::bary_mass_from_center_gm1(real_t gm1c) 
const -> real_t
{
  return mb_gm1(gm1c);
}

auto star_seq_impl::circ_radius_from_center_gm1(real_t gm1c) 
const -> real_t
{
  return rc_gm1(gm1c);
}

auto star_seq_impl::moment_inertia_from_center_gm1(real_t gm1c) 
const -> real_t
{
  return mi_gm1(gm1c);
}

auto star_seq_impl::lambda_tidal_from_center_gm1(real_t gm1c) 
const -> real_t
{
  return lt_gm1(gm1c);
}

auto star_seq_impl::range_center_gm1() const -> range_t
{
  return rg_gm1;
}

auto star_seq_impl::contains_gm1(real_t gm1c) const -> bool
{ 
  return range_center_gm1().contains(gm1c);
}



auto star_seq_impl::units_to_SI() const -> const units&
{
  return u;
}



star_branch_impl::star_branch_impl(range_t rg_gm1_, 
           spline_t xg_mg_, real_t gm1_ref_, bool incl_max_, units u_)
: rg_gm1{rg_gm1_}, xg_mg{std::move(xg_mg_)}, 
  gm1_ref{gm1_ref_}, incl_max{incl_max_}, u{u_}
{}


auto star_branch_impl::center_gm1_from_grav_mass(real_t mg) 
const -> real_t
{
  // limit above 0 or roundoff errors could cause NAN from sqrt
  real_t xg{ std::max(0., xg_mg(mg)) };
  return gm1_from_xg(xg);
}


auto star_branch_impl::range_grav_mass() const -> range_t
{
  return xg_mg.range_x();  
}

auto star_branch_impl::contains_grav_mass(real_t mg) const -> bool
{ 
  return range_grav_mass().contains(mg);
}

auto star_branch_impl::range_center_gm1() const -> range_t
{
  return rg_gm1;  
}

auto star_branch_impl::contains_gm1(real_t gm1c) const -> bool
{ 
  return range_center_gm1().contains(gm1c);
}

auto star_branch_impl::includes_maximum() const -> bool
{
  return incl_max;
}

auto star_branch_impl::grav_mass_maximum() const -> real_t
{
  return range_grav_mass().max();
}

//~ auto star_branch_impl::bary_mass_maximum() const -> real_t
//~ {
  //~ return bary_mass_from_center_gm1(center_gm1_maximum());  
//~ }

auto star_branch_impl::center_gm1_maximum() const -> real_t
{
  return rg_gm1.max();
}

auto star_branch_impl::xg_from_gm1(real_t gm1, real_t gm1_ref) 
-> real_t
{
  return std::pow(gm1_ref - gm1, 2);
}
  
auto star_branch_impl::gm1_from_xg(real_t xg) const -> real_t
{
  return gm1_ref - std::sqrt(xg);
}

}  // namespace detail




auto star_seq::valid() const -> const impl_t&
{
  assert(pimpl);
  return *pimpl;
}

auto star_seq::grav_mass_from_center_gm1(real_t gm1c) const 
-> real_t
{
  auto v = valid();
  return v.contains_gm1(gm1c) ? v.grav_mass_from_center_gm1(gm1c)
                              : numeric_limits<real_t>::quiet_NaN();
}

auto star_seq::bary_mass_from_center_gm1(real_t gm1c) const 
-> real_t
{
  auto v = valid();
  return v.contains_gm1(gm1c) ? v.bary_mass_from_center_gm1(gm1c)
                              : numeric_limits<real_t>::quiet_NaN();
  ;
}

auto star_seq::circ_radius_from_center_gm1(real_t gm1c) const 
-> real_t
{
  auto v = valid();
  return v.contains_gm1(gm1c) ? v.circ_radius_from_center_gm1(gm1c)
                              : numeric_limits<real_t>::quiet_NaN();
}

auto star_seq::moment_inertia_from_center_gm1(real_t gm1c) const 
-> real_t
{
  auto v = valid();
  return v.contains_gm1(gm1c) ? v.moment_inertia_from_center_gm1(gm1c)
                              : numeric_limits<real_t>::quiet_NaN();
}

auto star_seq::lambda_tidal_from_center_gm1(real_t gm1c) const 
-> real_t
{
  auto v = valid();
  return v.contains_gm1(gm1c) ? v.lambda_tidal_from_center_gm1(gm1c)
                              : numeric_limits<real_t>::quiet_NaN();
}

auto star_seq::range_center_gm1() const -> range_t
{
  return valid().range_center_gm1();
}

auto star_seq::contains_gm1(real_t gm1c) const -> bool
{
  return range_center_gm1().contains(gm1c);
}

auto star_seq::units_to_SI() const -> const units&
{
  return valid().units_to_SI();
}


auto star_branch::valid() const -> const impl_t&
{
  assert(pimpl);
  return *pimpl;
}

auto star_branch::center_gm1_from_grav_mass(real_t mg) const 
-> real_t
{
  auto v = valid();
  if (!v.contains_grav_mass(mg)) 
  {
    return numeric_limits<real_t>::quiet_NaN(); 
  }
  real_t gm1{ v.center_gm1_from_grav_mass(mg) };
  //make sure to always stay in range despite roundoff errors
  return v.range_center_gm1().limit_to(gm1);
}

auto star_branch::bary_mass_from_grav_mass(real_t mg) const 
-> real_t
{
  real_t gm1{ center_gm1_from_grav_mass(mg) };
  return star_seq::bary_mass_from_center_gm1(gm1);
}

auto star_branch::circ_radius_from_grav_mass(real_t mg) const 
-> real_t
{
  real_t gm1{ center_gm1_from_grav_mass(mg) };
  return star_seq::circ_radius_from_center_gm1(gm1);
}

auto star_branch::moment_inertia_from_grav_mass(real_t mg) const 
-> real_t
{
  real_t gm1{ center_gm1_from_grav_mass(mg) };
  return star_seq::moment_inertia_from_center_gm1(gm1);
}

auto star_branch::lambda_tidal_from_grav_mass(real_t mg) const 
-> real_t
{
  real_t gm1{ center_gm1_from_grav_mass(mg) };
  return star_seq::lambda_tidal_from_center_gm1(gm1);
}


auto star_branch::range_center_gm1() const -> range_t
{
  return valid().range_center_gm1();
}

auto star_branch::contains_gm1(real_t gm1c) const -> bool
{
  return valid().contains_gm1(gm1c);
}

auto star_branch::range_grav_mass() const -> range_t
{
  return valid().range_grav_mass();
}

auto star_branch::contains_grav_mass(real_t mg) const -> bool
{
  return valid().contains_grav_mass(mg);
}

auto star_branch::includes_maximum() const -> bool
{
  return valid().includes_maximum();
}

auto star_branch::grav_mass_maximum() const -> real_t
{
  return valid().grav_mass_maximum();
}

auto star_branch::bary_mass_maximum() const -> real_t
{
  return star_seq::bary_mass_from_center_gm1(center_gm1_maximum());
}

auto star_branch::center_gm1_maximum() const -> real_t
{
  return valid().center_gm1_maximum();
}

auto star_branch::grav_mass_from_center_gm1(real_t gm1c) const 
-> real_t
{
  return contains_gm1(gm1c) 
             ? star_seq::grav_mass_from_center_gm1(gm1c)
             : numeric_limits<real_t>::quiet_NaN();
}

auto star_branch::bary_mass_from_center_gm1(real_t gm1c) 
const -> real_t
{
  return contains_gm1(gm1c) 
             ? star_seq::bary_mass_from_center_gm1(gm1c)
             : numeric_limits<real_t>::quiet_NaN();
}

auto star_branch::circ_radius_from_center_gm1(real_t gm1c) 
const -> real_t
{
  return contains_gm1(gm1c) 
             ? star_seq::circ_radius_from_center_gm1(gm1c)
             : numeric_limits<real_t>::quiet_NaN();
}

auto star_branch::moment_inertia_from_center_gm1(real_t gm1c) 
const -> real_t
{
  return contains_gm1(gm1c) 
             ? star_seq::moment_inertia_from_center_gm1(gm1c)
             : numeric_limits<real_t>::quiet_NaN();
}

auto star_branch::lambda_tidal_from_center_gm1(real_t gm1c) 
const -> real_t
{
  return contains_gm1(gm1c) 
             ? star_seq::lambda_tidal_from_center_gm1(gm1c)
             : numeric_limits<real_t>::quiet_NaN();
}

}

//~ namespace {


template <class F, class R>
auto bracket_maximum(const F& f, R guess, EOS_Toolkit::interval<R> bnd, 
                     int max_step=100, double search_fac=1.5) 
-> EOS_Toolkit::interval<R>
{
  assert(bnd.max() > 0);
  assert(guess > 0);
  assert(search_fac > 1.);

  R xg{ bnd.limit_to(guess) };
  R x2{ bnd.limit_to(xg * search_fac) };
  R x0{ bnd.limit_to(x2 / (search_fac*search_fac)) };
  R x1{ sqrt(x2 * x0) };
  
  
  auto f0{ f(x0) };
  auto f1{ f(x1) };
  auto f2{ f(x2) };
  
  while(--max_step > 0) {
    if ((f0 < f1) && (f1 > f2)) {
      return {x0, x2};
    }
    if (f2 >= f0) {
      x0 = x1; f0 = f1;
      x1 = x2; f1 = f2;
      x2 = x2 * search_fac;
      if (x2 >= bnd.max()) return {x0, bnd.max()};
      f2 = f(x2);
    } 
    else {
      x2 = x1; f2 = f1;
      x1 = x0; f1 = f0;
      x0 = x0 / search_fac;
      if (x0 <= bnd.min()) return {bnd.min(), x2};
      f0 = f(x0);
    }
  }
  throw std::runtime_error("Maximum search failed (too many steps)");
}


template <class F, class R>
auto bracket_root(const F& f, R guess, EOS_Toolkit::interval<R> bnd, 
                  unsigned int max_step=100, double search_fac=1.5) 
-> EOS_Toolkit::interval<R>
{
  
  assert( search_fac > 1 );
  assert( bnd.max() > 0 );
  assert( guess > 0 );
  
  
  R x1{ bnd.limit_to(guess) };
  R x0{ bnd.limit_to(x1 / search_fac)};
  
  auto f0{ f(x0) };
  auto f1{ f(x1) };
  
  while(--max_step > 0) {
    if  (((f0 < 0) && (0 < f1)) || ((f1 < 0) && (0 < f0))) {
      return {x0, x1};
    }
    if (((f0 >= 0) && (f1 > f0)) || ( (f0 <= 0) && (f1 < f0) )) {
      x0 = x0 / search_fac;
      if (x0 <= bnd.min()) {
        if (f(bnd.min())*f1 <= 0) return {bnd.min(),x1};
        throw std::runtime_error("Root bracket failed (out of bounds)");
      }
      f0 = f(x0);
    } 
    else {
      x1 = x1 * search_fac;
      if (x1 >= bnd.max()) {
        if (f(bnd.max())*f0 <= 0) return {x0,bnd.max()};
        throw std::runtime_error("Root bracket failed (out of bounds)");
      }
      f1 = f(x1);
    }
  }
  throw std::runtime_error("Root bracket failed (too many steps)");
}

template <class F, class R, class V>
auto bracket_value(const F& f, V v, R guess, EOS_Toolkit::interval<R> bnd, 
                   unsigned int max_step=100, double search_fac=1.5) 
-> EOS_Toolkit::interval<R>
{
  
  auto f2 = [&] (R x) {
    return f(x) - v;
  };
  
  return bracket_root(f2, guess, bnd, max_step, search_fac);
}

template <class F, class R>
auto find_maximum(const F& f, EOS_Toolkit::interval<R> bnd, const int bits=40, 
                  unsigned int max_steps=100) 
-> R
{
  boost::uintmax_t iters{ max_steps };
  const auto r{ 
    boost::math::tools::brent_find_minima(
                               [&] (EOS_Toolkit::real_t x) {return -f(x);}, 
                               bnd.min(), bnd.max(), 
                               bits, iters)
  };
  
  if (iters >= max_steps)
  {
    throw std::runtime_error("maximum not found");
  }
  
  return r.first;
}



namespace EOS_Toolkit {

auto make_tov_seq_impl(eos_barotr eos, tov_acc_simple acc, 
                       interval<real_t> rg_gm1, unsigned int num_samp)
-> std::shared_ptr<detail::star_seq_impl>
{
  assert(num_samp>5);
  
  std::vector<real_t> mg(num_samp), mb(num_samp), rc(num_samp), 
                      mi(num_samp), lt(num_samp);
  
  for (unsigned int i=0; i < num_samp; ++i)
  {
    const real_t w{ i/double(num_samp-1) };
    const real_t gm1{ 
      eos.range_gm1().limit_to( rg_gm1.min() + w * rg_gm1.length() ) 
    };
    const real_t rhoc{ eos.at_gm1(gm1).rho() };
    auto tov{ get_tov_star_properties(eos, rhoc, acc) };
    
    mg[i] = tov.grav_mass();
    mb[i] = tov.bary_mass();
    rc[i] = tov.circ_radius();
    mi[i] = tov.moment_inertia();
    lt[i] = tov.deformability().lambda;
  }
  
  return detail::star_seq_impl::from_vector(std::move(mg),
                                 std::move(mb), std::move(rc),
                                 std::move(mi), std::move(lt),
                                 rg_gm1, eos.units_to_SI());

}

auto make_tov_seq(eos_barotr eos, tov_acc_simple acc, 
                       interval<real_t> rg_gm1, unsigned int num_samp)
-> star_seq
{
  return star_seq(make_tov_seq_impl(eos, acc, rg_gm1, num_samp));
}



auto make_tov_branch_impl(const detail::star_seq_impl& seq,
                          interval<real_t> rg_gm1, real_t gm1_ref,
                          std::size_t nsamp,
                          std::size_t tmp_nsamp, bool incl_max)
-> std::shared_ptr<detail::star_branch_impl>
{
  std::vector<real_t> tmp_mg(tmp_nsamp), tmp_xg(tmp_nsamp);
  
  for (unsigned int i=0; i < tmp_nsamp; ++i)
  {
    const real_t w{ i/double(tmp_nsamp-1) };
    const real_t gm1{ rg_gm1.min() + w * rg_gm1.length() };
    tmp_xg[i] = detail::star_branch_impl::xg_from_gm1(gm1, gm1_ref);
    tmp_mg[i] = seq.grav_mass_from_center_gm1(gm1);
  }
  
  auto tmp_xg_mg = make_interpol_pchip_spline(tmp_mg, tmp_xg);
  
  interval<real_t> rg_mg{
    seq.grav_mass_from_center_gm1(rg_gm1.min()),
    seq.grav_mass_from_center_gm1(rg_gm1.max()) 
  };
  
  std::vector<real_t> xg(nsamp);
  for (unsigned int i=0; i < nsamp; ++i)
  {
    const real_t w{ i/double(nsamp-1) };
    const real_t mg{ rg_mg.min() + w * rg_mg.length() };
    xg[i] = tmp_xg_mg(mg);
  }
  auto xg_mg = make_interpol_regspl(xg, rg_mg);

  return make_shared<detail::star_branch_impl>(rg_gm1, xg_mg, 
                                       gm1_ref, incl_max, 
                                       seq.units_to_SI());
}


auto make_tov_branch_stable(eos_barotr eos, tov_acc_simple acc,
                          real_t mgrav_min,
                          unsigned int num_samp,
                          real_t gm1_initial, real_t max_margin)
-> star_branch
{
  const unsigned int oversamp=2;
  const unsigned int tmp_subsamp=10;
  
  if (max_margin <= 0) {
    throw std::invalid_argument("Margin for true maximum must be" 
                                "positive");
  }
  
  auto f = [&] (real_t gm1c) {
    real_t rhoc{ eos.range_rho().limit_to(eos.at_gm1(gm1c).rho()) };
    auto tov{ get_tov_star_properties(eos, rhoc, acc) };
    return tov.grav_mass();
  };
  auto bracket_max{ bracket_maximum(f, gm1_initial, eos.range_gm1()) };
  
  auto bracket_low{ 
    bracket_value(f, mgrav_min, 
                  std::min(gm1_initial, bracket_max.min()),
                  eos.range_gm1()) 
  };
  
  interval<real_t> rg_seq_gm1{bracket_low.min(), bracket_max.max()};
  
  auto sqimpl{ make_tov_seq_impl(eos, acc, rg_seq_gm1, num_samp) };
  
  auto f2 = [&] (real_t gm1c) {
    return sqimpl->grav_mass_from_center_gm1(gm1c);
  };
  
  const real_t gm1_max{ find_maximum(f2, bracket_max) };
    
  const bool incl_max{ 
    eos.range_gm1().contains(gm1_max * (1.0 + max_margin))
  };
  
  
  interval<real_t> rg_gm1{sqimpl->range_center_gm1().min(), 
                          std::min(sqimpl->range_center_gm1().max(), 
                                   gm1_max)};

  auto bimpl{ 
    make_tov_branch_impl(*sqimpl, rg_gm1, gm1_max,
                num_samp*oversamp, num_samp*tmp_subsamp, incl_max) 
  };
  
  return star_branch{sqimpl, bimpl};
}


auto make_star_seq(std::vector<real_t> mg, 
       std::vector<real_t> mb, std::vector<real_t> rc, 
       std::vector<real_t> mi, std::vector<real_t> lt, 
       star_seq::range_t rg_gm1, units u)
-> star_seq
{
  return star_seq{
    detail::star_seq_impl::from_vector(std::move(mg),
                                 std::move(mb), std::move(rc),
                                 std::move(mi), std::move(lt),
                                 rg_gm1, u)
  };
}



}


