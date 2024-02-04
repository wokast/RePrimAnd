#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>
#include <stdexcept>
#include <deque>
#include <cmath>
#include <memory>
#include <algorithm>
#include <limits>
#include <iostream>
#include <iomanip>
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
  
  const auto accs{ star_acc_simple(false, false, acc) };
  
  auto fmin = [&] (real_t rho_c) {
    auto tov{ 
      get_tov_properties(eos, rho_c, accs) 
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
  const auto accs{ star_acc_simple(false, false, acc) };
  
  auto froot = [&] (real_t rho_c) {
    auto tov{ get_tov_properties(eos, rho_c, accs) };
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
  lt_gm1{std::move(lt_gm1_)}, rg_gm1{rg_gm1_}, 
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
  if (!mg_gm1.range_x().contains(rg_gm1))
  {
    throw std::runtime_error("Attempt to create star sequence with "
                      "insufficient range for M_g(g-1)");
  } 
  if (!mb_gm1.range_x().contains(rg_gm1))
  {
    throw std::runtime_error("Attempt to create star sequence with "
                      "insufficient range for M_b(g-1)");
  } 
  if (!rc_gm1.range_x().contains(rg_gm1))
  {
    throw std::runtime_error("Attempt to create star sequence with "
                      "insufficient range for R(g-1)");
  } 
  if (!mi_gm1.range_x().contains(rg_gm1))
  {
    throw std::runtime_error("Attempt to create star sequence with "
                      "insufficient range for I(g-1)");
  } 
  if (!lt_gm1.range_x().contains(rg_gm1))
  {
    throw std::runtime_error("Attempt to create star sequence with "
                      "insufficient range for Lambda(g-1)");
  } 
}


star_seq_impl::star_seq_impl(spline_t mg_gm1_, spline_t mb_gm1_, 
                spline_t rc_gm1_, spline_t mi_gm1_, 
                spline_t lt_gm1_, units u_)
: star_seq_impl(mg_gm1_, mb_gm1_, rc_gm1_, mi_gm1_, lt_gm1_, 
                intersect(mg_gm1_.range_x(), mb_gm1_.range_x(), 
                          rc_gm1_.range_x(), mi_gm1_.range_x(), 
                          lt_gm1_.range_x()), 
                u_)
{}

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


auto star_seq_impl::from_logspaced_samples(std::vector<real_t> mg, 
       std::vector<real_t> mb, std::vector<real_t> rc, 
       std::vector<real_t> mi, std::vector<real_t> lt, 
       range_t rg_gm1, units u)
-> std::shared_ptr<star_seq_impl>
{

  auto mg_gm1 = make_interpol_logspl(std::move(mg), rg_gm1);
  auto mb_gm1 = make_interpol_logspl(std::move(mb), rg_gm1);
  auto rc_gm1 = make_interpol_logspl(std::move(rc), rg_gm1);
  auto mi_gm1 = make_interpol_logspl(std::move(mi), rg_gm1);
  auto lt_gm1 = make_interpol_logspl(std::move(lt), rg_gm1);
  
  
  return make_shared<star_seq_impl>(mg_gm1, mb_gm1, rc_gm1,
                                    mi_gm1, lt_gm1, u);
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




namespace EOS_Toolkit {

auto make_star_seq_impl(
       std::function<spherical_star_properties(real_t)> solver,
       interval<real_t> rg_gm1, units gu, unsigned int num_samp)
-> std::shared_ptr<detail::star_seq_impl>
{
  assert(num_samp>5);
  
  std::vector<real_t> mg(num_samp), mb(num_samp), rc(num_samp), 
                      mi(num_samp), lt(num_samp);
  
  auto rg_lggm1{ log(rg_gm1) };
  for (unsigned int i=0; i < num_samp; ++i)
  {
    const real_t w{ i/double(num_samp-1) };
    const real_t lggm1{ 
      rg_lggm1.limit_to( rg_lggm1.min() + w * rg_lggm1.length() ) 
    };
    const real_t gm1{ std::exp(lggm1) };
    auto tov{ solver(gm1) };
    
    mg[i] = tov.grav_mass();
    mb[i] = tov.bary_mass();
    rc[i] = tov.circ_radius();
    mi[i] = tov.moment_inertia();
    lt[i] = tov.deformability().lambda;
  }
  
  return detail::star_seq_impl::from_logspaced_samples(std::move(mg),
                                 std::move(mb), std::move(rc),
                                 std::move(mi), std::move(lt),
                                 rg_gm1, gu);
}

auto make_tov_seq_impl(eos_barotr eos, interval<real_t> rg_gm1, 
                       const star_accuracy_spec acc,
                       unsigned int num_samp)
-> std::shared_ptr<detail::star_seq_impl>
{
  if (!eos.range_gm1().contains(rg_gm1)) 
  {
    throw std::runtime_error("make_tov_seq_impl: requested range "
                             "exceeds EOS validity range");
  }
  auto solver = [&] (real_t gm1) {
    const real_t rhoc{ eos.at_gm1(gm1).rho() };
    return  get_tov_properties(eos, rhoc, acc);
  };
  
  return make_star_seq_impl(solver, rg_gm1, eos.units_to_SI(), 
                            num_samp);
}


auto make_star_seq(
      std::function<spherical_star_properties(real_t)> solver,
      interval<real_t> rg_gm1, units u, unsigned int num_samp)
-> star_seq
{
  return star_seq(make_star_seq_impl(solver, rg_gm1, u, num_samp));
}


auto make_tov_seq(eos_barotr eos, interval<real_t> rg_gm1, 
                  const star_accuracy_spec acc, unsigned int num_samp)
-> star_seq
{
  return star_seq(make_tov_seq_impl(eos, rg_gm1, acc,  num_samp));
}


auto make_tov_branch_impl(const detail::star_seq_impl& seq,
                          interval<real_t> rg_gm1, real_t gm1_ref,
                          std::size_t tmp_nsamp, bool incl_max)
-> std::shared_ptr<detail::star_branch_impl>
{
  assert(gm1_ref>=rg_gm1);
  assert(0. < rg_gm1);
  
  std::vector<real_t> tmp_mg, tmp_xg;
  tmp_mg.reserve(tmp_nsamp);
  tmp_xg.reserve(tmp_nsamp);
  auto rg_lggm1{ log(rg_gm1) };
  
  for (unsigned int i=0; i < tmp_nsamp-1; ++i)
  {
    const real_t w{ i/double(tmp_nsamp-1) };
    const real_t lggm1{ rg_lggm1.min() + w * rg_lggm1.length() };
    const real_t gm1{ std::exp(lggm1) };
    tmp_xg.push_back( detail::star_branch_impl::xg_from_gm1(gm1, gm1_ref) );
    tmp_mg.push_back( seq.grav_mass_from_center_gm1(gm1) );
  }
  tmp_xg.push_back( detail::star_branch_impl::xg_from_gm1(rg_gm1.max(), gm1_ref) );
  tmp_mg.push_back( seq.grav_mass_from_center_gm1(rg_gm1.max()) );
  
  auto tmp_xg_mg = make_interpol_pchip_spline(tmp_mg, tmp_xg);

  return make_shared<detail::star_branch_impl>(rg_gm1, tmp_xg_mg, 
                                       gm1_ref, incl_max, seq.units_to_SI());
}

namespace {

class increment_limit {
   std::size_t count;
   const std::string msg;
   
   public:
   increment_limit(std::size_t count_, std::string msg_) 
   : count{count_}, msg{msg_} {}
   increment_limit& operator++(int) {
     if (--count == 0) throw std::runtime_error(msg);
     return *this;
   }
};  

template<class F, class I, class T=decltype(std::declval<F>()(I()))>
auto filter_indices_ascending(F f, const I imin, const I imax, const I iend)
-> std::pair< std::vector<I>, std::vector<I> >
{
  const T fmin{ f(imin) };
  const T fmax{ f(imax) };
  assert(fmin < fmax);
  assert(imin>=0);
  assert(imax>=imin);
  assert(imax<=iend);
  
  std::pair <std::vector<I>, std::vector<I> > k;
  
  if (imin>0) 
  {  
    k.first.push_back(0);
  }
  k.first.push_back(imin);
  for (auto i=imin; i < imax; ++i)
  {
    const T fi{ f(i) };
    if ((f(k.first.back()) < fi) && (fi < fmax)) 
    {
      k.first.push_back( i );
    }
  }
  k.first.push_back(imax);
  if (imax<iend) 
  {  
    k.first.push_back(iend);
  }
  
  std::vector<I> k_lo_rev;
  
  
  k_lo_rev.push_back(imax);
  for (auto i=imax; i > imin; --i)
  {
    const T fi{ f(i) };
    if ((f(k_lo_rev.back()) > fi)  && (fi > fmin))
    {
      k_lo_rev.push_back( i );
    }
  }
  k_lo_rev.push_back(imin);
  
  if (imin>0) 
  {  
    k.second.push_back(0);
  }
  for (auto i=k_lo_rev.rbegin(); i<k_lo_rev.rend(); ++i) 
  {
    k.second.push_back(*i);
  }
  if (imax<iend) 
  {  
    k.second.push_back(iend);
  }
  return k;
}

template<class X, class Y, class I, class T=decltype(std::declval<Y>()(I()))>
auto interp_seq (X x, Y y, std::vector<I> k)
-> interpolator
{
  std::vector<T> vx,vy;
  for (auto i : k)
  { 
    vx.push_back(x(i));
    vy.push_back(y(i));
  }
  return make_interpol_pchip_spline(vx, vy);
}

template<class F, class G, class A=real_t, class V=A>
auto func_average(const F& f, const G& g) 
-> std::function<V(A)>
{
  return [&] (A x) {return (f(x) + g(x)) / 2.;};
}

template<class X, class Y, class I>
auto make_special_monot_interp(const X& xi, const Y& yi, 
               const interval<real_t> rg_x_mono, 
               const interval<real_t> rg_x_full, 
               const std::vector<I>& idx_hi, const std::vector<I>& idx_lo,
               const std::size_t nresamp) 
-> interpolator
{
  
  auto rg_lgx_mono{ log(rg_x_mono) };
  auto rg_lgx_full{ log(rg_x_full) };
  const real_t dlgx{ rg_lgx_mono.length() / (nresamp-1) };
  
  const int nrbnd{ 
    static_cast<int>(floor((rg_lgx_full.max() - rg_lgx_mono.max()) / dlgx)) 
  };
  const int nlbnd{ 
    static_cast<int>(floor((rg_lgx_mono.min() - rg_lgx_full.min()) / dlgx)) 
  };
  
  assert(nlbnd>=0);
  assert(nrbnd>=0);
  
  interval<real_t> rg_lgx_adj{ 
    rg_lgx_mono.min() - nlbnd*dlgx, 
    rg_lgx_mono.max() + nrbnd*dlgx
  };
  
  auto hi{ interp_seq(xi, yi, idx_hi) };
  auto lo{ interp_seq(xi, yi, idx_lo) };
  auto av{ func_average(hi, lo) };
  
  return make_interpol_logspl(av, exp(rg_lgx_adj), nresamp+nlbnd+nrbnd);
}


struct sampled_stable_branch { 
  std::deque<spherical_star_properties> seq;
  std::ptrdiff_t imin, imax, ilast;
};

auto scout_stable_branch(
  std::function<spherical_star_properties(real_t)> solver,
  const interval<real_t> val_rg_gm1, 
  const real_t acc_mg,
  const real_t mg_cut_low_rel,
  const real_t mg_cut_low_abs,
  const real_t gm1_init, 
  const real_t gm1_step) 
-> sampled_stable_branch
{
  const real_t tolf{ 1. + acc_mg }; 
  const real_t gm1fac{ 1. + gm1_step };

  if (gm1_init <= 0.0) {
    throw std::invalid_argument("Star branch search: start enthalpy " 
                                "violates g-1 > 0");
  }
  const real_t gm1_start{ val_rg_gm1.limit_to(gm1_init) };
  
  if (gm1fac <= 1.0) {
    throw std::invalid_argument("Star branch search: enthalpy step size" 
                                "must be strictly positive");
  }
  if (mg_cut_low_rel >= 1) {
    throw std::invalid_argument("Star branch search: mass cutoff " 
                                "factor must be < 1");
  }
  
  if (tolf <= 1) {
    throw std::invalid_argument("Star branch search: allowed TOV solver " 
                                "relative error must be strictly positive");
  }

  std::deque<spherical_star_properties> sq;
  
  
  const auto tov_initial{ solver(gm1_start) };
  const real_t mg_initial{ tov_initial.grav_mass() };
  
  sq.push_back (tov_initial);
  real_t mg_max{ mg_initial };
  
  {
    real_t gm1_raise{ gm1_start * gm1fac };
    while ((sq.back().grav_mass() * tolf >= mg_max / tolf ) 
             && val_rg_gm1.contains(gm1_raise))  
    {
      sq.push_back(solver(gm1_raise));
      mg_max = std::max(mg_max, sq.back().grav_mass());
      gm1_raise *= gm1fac;
    }
    
    for (int i=0; i<3; ++i) 
    {
      if (!val_rg_gm1.contains(gm1_raise)) break;  
      sq.push_back(solver(gm1_raise));
      gm1_raise *= gm1fac;
    }
  
  }

  real_t gm1_fall{ gm1_start / gm1fac };
  while ((sq.front().grav_mass() * tolf >= mg_max / tolf ) 
           && val_rg_gm1.contains(gm1_fall))
  {
    sq.push_front(solver(gm1_fall));
    mg_max = std::max(mg_max, sq.front().grav_mass());
    gm1_fall /= gm1fac;
  } 

  if (mg_cut_low_abs >= mg_max) {
    throw std::invalid_argument("Star branch search: lower mass cutoff " 
                                "is above maxmimum mass");
  }
  
  const real_t mg_cut_min{ std::max(mg_cut_low_abs, mg_max * mg_cut_low_rel) };
  
  real_t mg_min{ sq.front().grav_mass() };
  while ((sq.front().grav_mass() / tolf <= mg_min * tolf ) 
           && (sq.front().grav_mass() * tolf >= mg_cut_min) 
           && val_rg_gm1.contains(gm1_fall)) 
  {
    sq.push_front(solver(gm1_fall));
    mg_min = std::min(mg_min, sq.front().grav_mass());
    gm1_fall /= gm1fac;
  } 
  
  for (int i=0; i<3; ++i) 
  {
    if (!val_rg_gm1.contains(gm1_fall)) break;  
    sq.push_front(solver(gm1_fall));
    gm1_fall /= gm1fac;
  }
  
  
  auto itmin = sq.end();
  auto itmax = sq.end();
  for (auto i=sq.begin(); i!=sq.end(); ++i) {
    if (i->grav_mass() == mg_max) itmax=i;
    if (i->grav_mass() == mg_min) itmin=i;
  }
  
  assert(sq.size() > 2);
  assert(itmin != sq.end());
  assert(itmax != sq.end());
  assert(itmin != itmax);
  
  std::ptrdiff_t imin{ itmin - sq.begin() };
  std::ptrdiff_t imax{ itmax - sq.begin() };
  std::ptrdiff_t ilast{ std::prev(sq.end()) - sq.begin() };
  
  const std::ptrdiff_t nsamp_min{ 10 };
  
  if (imax  < imin + nsamp_min)  
  {
    const real_t gm1_init_new{ 
      (sq[imin].center_gm1() + sq[imax].center_gm1()) / 2. 
    };
    const real_t gm1_step_new{ gm1_step / nsamp_min };
    return scout_stable_branch(solver, val_rg_gm1, acc_mg, 
                 mg_cut_low_rel, mg_cut_low_abs,
                 gm1_init_new, gm1_step_new);  
    
  }
  
  return sampled_stable_branch{sq, imin, imax, ilast};
}




auto improve_tovseq_extremum(
  const spherical_star_properties& s0,
  const spherical_star_properties& s1,
  const spherical_star_properties& s2,
  std::function<spherical_star_properties(real_t)> solver
)   
-> spherical_star_properties
{
  const real_t y0{ s0.grav_mass() };
  const real_t y1{ s1.grav_mass() };
  const real_t y2{ s2.grav_mass() };
  
  const bool ismax{ (y0<=y1) && (y2<=y1) && ((y0<y1) || (y2<y1)) };
  const bool ismin{ (y0>=y1) && (y2>=y1) && ((y0>y1) || (y2>y1)) };
  
  if (!(ismin || ismax)) return s1;
  
  const real_t x0{ s0.center_gm1() };
  const real_t x1{ s1.center_gm1() };
  const real_t x2{ s2.center_gm1() };
  const real_t d01{ x0 - x1 };
  const real_t d02{ x0 - x2 };
  const real_t d12{ x1 - x2 };
  
  const real_t z0{ y0 / (d01*d02) };
  const real_t z1{ -y1 / (d01*d12) };
  const real_t z2{ y2 / (d02*d12) };
  
  const real_t xmax{ (z0*(x1+x2) + z1*(x0+x2) + z2*(x1+x0)) / (2.*(z0+z1+z2)) };
  assert(xmax>x0);
  assert(xmax<x2);
  
  auto sa{ solver(xmax) };
  
  assert(sa.center_gm1() > x0);
  assert(sa.center_gm1() < x2);
  
  if ((ismax && (sa.grav_mass() <= y1)) || (ismin && (sa.grav_mass() >= y1)))
  {
    return s1;
  }
  return sa;
}

} // anonymous namespace


auto make_star_branch_stable(
        std::function<spherical_star_properties(real_t)> solver,
        const interval<real_t> val_rg_gm1, const real_t acc_mg,
        const units gu, const real_t mg_cut_low_rel,
        const real_t mg_cut_low_abs,
        const real_t gm1_initial, const real_t gm1_step, 
        const real_t max_margin)
-> star_branch
{
  const unsigned int tmp_subsamp{ 5 };
  const std::ptrdiff_t min_nresamp{ 10 };
  
  auto scb{ scout_stable_branch(solver, val_rg_gm1, acc_mg,
                        mg_cut_low_rel, mg_cut_low_abs, gm1_initial, gm1_step) };
    

  const std::ptrdiff_t imin{ scb.imin };
  const std::ptrdiff_t imax{ scb.imax };
  const std::ptrdiff_t ilast{ scb.ilast };

  

  if ((imax > 0) && (imax < ilast))
  {
    scb.seq[imax] = improve_tovseq_extremum(scb.seq[imax-1], 
                                            scb.seq[imax], 
                                            scb.seq[imax+1], solver);
  }
  if ((imin > 0) && (imin < ilast))
  {
    scb.seq[imin] = improve_tovseq_extremum(scb.seq[imin-1], 
                                            scb.seq[imin], 
                                            scb.seq[imin+1], solver);
  }


  auto qati = [&scb] (real_t (spherical_star_properties::*fun)() const) {
    return [fun, &scb] (std::size_t i) {
      spherical_star_properties t{ scb.seq.at(i) };
      return (t.*fun)();
    };
  };

  auto mgati = qati(&spherical_star_properties::grav_mass);
  
  auto kasc{ filter_indices_ascending(mgati, imin, imax, ilast) };
  
  auto gm1ati = qati(&spherical_star_properties::center_gm1);
  
  auto ltati = [&scb] (std::size_t i) {
    return scb.seq.at(i).deformability().lambda;
  };


  interval<real_t> rg_gm1_inc{gm1ati(imin), gm1ati(imax)};
  interval<real_t> rg_gm1_full{gm1ati(0), gm1ati(ilast)};
  
  std::ptrdiff_t num_samp{ std::max(min_nresamp, imax - imin + 1)};
  
  auto mkip = [&] (std::function<real_t(std::size_t)> yi) 
  {
    return make_special_monot_interp(gm1ati, yi, rg_gm1_inc, rg_gm1_full, 
                                     kasc.first, kasc.second, num_samp); 
  };
  
  auto mggm1 = mkip(mgati);  
  auto mbgm1 = mkip(qati(&spherical_star_properties::bary_mass));
  auto rcgm1 = mkip(qati(&spherical_star_properties::circ_radius));
  auto migm1 = mkip(qati(&spherical_star_properties::moment_inertia));
  auto ltgm1 = mkip(ltati);
  
  auto sqimpl{
    std::make_shared<detail::star_seq_impl>(mggm1, mbgm1, rcgm1, migm1, ltgm1, 
                                            mggm1.range_x(), gu)
  };
      
  
  const bool incl_max{ 
    val_rg_gm1.contains(rg_gm1_inc.max() * (1.0 + max_margin))
  };
  
  real_t gm1_ref{ 
    rg_gm1_inc.max() * (incl_max ? 1 : 2)
  };

  
  auto bimpl{ 
    make_tov_branch_impl(*sqimpl, rg_gm1_inc, gm1_ref,
                num_samp*tmp_subsamp, incl_max) 
  };
  
  return star_branch{sqimpl, bimpl};  
  
}


auto make_tov_branch_stable(eos_barotr eos, 
            const star_accuracy_spec acc,
            real_t mg_cut_low_rel, real_t mg_cut_low_abs, 
            real_t gm1_initial, 
            real_t gm1_step, real_t max_margin)
-> star_branch
{
  increment_limit fin(1000000, "TOV branch search seems stuck, aborting");
  
  auto solver = [&] (real_t gm1) {
    fin++;
    const real_t rhoc{ eos.at_gm1(gm1).rho() };
    auto tov = get_tov_properties(eos, rhoc, acc);
    assert(std::isfinite(tov.grav_mass()));
    assert(tov.grav_mass()>0);
    return tov;
  };
  
  return make_star_branch_stable(
      solver, eos.range_gm1(), 2*acc.acc_mass,
      eos.units_to_SI(), mg_cut_low_rel, mg_cut_low_abs,
      gm1_initial, gm1_step, max_margin
  );
  
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


