#include "eos_barotr_spline.h"
#include "eos_barotr_spline_impl.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <iomanip>


using namespace std;
using namespace EOS_Toolkit;
using namespace EOS_Toolkit::implementations;

using func_t  = std::function<real_t(real_t)>;

auto eos_barotr_spline::get_rggm1(const lgspl_t& eps_gm1, 
  const lglgspl_t& p_gm1, const lgspl_t& hm1_gm1, 
  const lglgspl_t& rho_gm1, const opt_t& temp_gm1, 
  const opt_t& efrac_gm1)
-> range
{
  range rg{ 
    intersect(eps_gm1.range_x(), p_gm1.range_x(), hm1_gm1.range_x(),
              rho_gm1.range_x())
  };
  
  if (temp_gm1) 
  {
    rg = intersect(rg, temp_gm1->range_x()); 
  }
  
  if (efrac_gm1) 
  {
    rg = intersect(rg, efrac_gm1->range_x()); 
  }

  return {0, rg.max()};
}

eos_barotr_spline::eos_barotr_spline(
    lglgspl_t gm1_rho_, lglgspl_t rho_gm1_,
    lgspl_t eps_gm1_, lglgspl_t p_gm1_, 
    lgspl_t hm1_gm1_, lgspl_t csnd_rho_, 
    opt_t temp_gm1_, opt_t efrac_gm1_, 
    bool isentropic_,
    const eos_barotr_gpoly& poly_)
: eos_barotr_impl{poly_.units_to_SI()},
  gm1_rho{ std::move(gm1_rho_) }, 
  eps_gm1{ std::move(eps_gm1_) },
  p_gm1{ std::move(p_gm1_) }, 
  hm1_gm1{ std::move(hm1_gm1_) },
  rho_gm1{ std::move(rho_gm1_) }, 
  csnd_rho{ std::move(csnd_rho_) }, 
  temp_gm1{ std::move(temp_gm1_) }, 
  efrac_gm1{ std::move(efrac_gm1_) }, 
  poly{ poly_ },  
  rggm1{ get_rggm1(eps_gm1, p_gm1, hm1_gm1, rho_gm1,
                   temp_gm1, efrac_gm1) },
  rgrho{ 0, rho_gm1(rggm1.max()) },
  gm1_low{ poly.range_gm1().max() },
  rho_low{ poly.range_rho().max() },
  min_h{ 1.0 + std::min(poly.hm1(0.), hm1_gm1.range_y().min()) },
  isentropic{isentropic_}
{ 
  
  if (!(eps_gm1.contains(gm1_low) &&
        p_gm1.contains(gm1_low) &&
        rho_gm1.contains(gm1_low) &&
        hm1_gm1.contains(gm1_low) &&
        (!temp_gm1 || temp_gm1->contains(gm1_low)) &&
        (!efrac_gm1 || efrac_gm1->contains(gm1_low)) ))
  {
    throw runtime_error("eos_barotr_spline: matching polytrope "
                        "outside sampled range for g-1"); 
  }
  
  
  if (!(gm1_rho.contains(rho_low) && csnd_rho.contains(rho_low)))
  {
    throw runtime_error("eos_barotr_spline: matching polytrope "
                        "outside sampled range for rho");     
  }
  
  if (rho_gm1.range_y().min() < 0.0) {
    throw runtime_error("eos_barotr_spline: negative mass density "
                        "in rho(gm1)");
  }
  if (csnd_rho.range_y().max() >= 1.0) {
    throw runtime_error("eos_barotr_spline: sound speed >= 1");
  }
  if (csnd_rho.range_y().min() < 0.0) {
    throw runtime_error("eos_barotr_spline: sound speed < 0");
  }
  if (p_gm1.range_y().min() < 0.0) {
    throw runtime_error("eos_barotr_spline: negative pressure");
  }
  if (gm1_rho.range_y().min() < 0.0) {
    throw runtime_error("eos_barotr_spline: encountered g < 1");
  }
  
  if (temp_gm1) {
    temp0    = (*temp_gm1)(rggm1.min());
    if (temp_gm1->range_y().min() < 0.0) {
      throw runtime_error("eos_barotr_spline: encountered negative "
                          "temperature");
    }
    zerotemp = temp_gm1->range_y().max() == 0;
  }
  
  if (zerotemp && !isentropic) {
    throw runtime_error("eos_barotr_spline: zero-temperature EOS must "
                        "be isentropic");
  }
  
  if (efrac_gm1) {
    efrac0    = (*efrac_gm1)(rggm1.min());
  }
}

real_t eos_barotr_spline::gm1_from_rho(real_t rho) const
{
  return (rho >= rho_low) ? gm1_rho(rho) : poly.gm1_from_rho(rho); 
}


real_t eos_barotr_spline::eps(real_t gm1) const
{
  return (gm1 >= gm1_low) ? eps_gm1(gm1) : poly.eps(gm1);
}


real_t eos_barotr_spline::press(real_t gm1) const
{
  return (gm1 >= gm1_low) ? p_gm1(gm1) : poly.press(gm1);
}

real_t eos_barotr_spline::rho(real_t gm1) const
{
  return (gm1 >= gm1_low) ? rho_gm1(gm1) : poly.rho(gm1);
}

real_t eos_barotr_spline::hm1(real_t gm1) const
{
  return (gm1 >= gm1_low) ? hm1_gm1(gm1)  : poly.hm1(gm1);
}

real_t eos_barotr_spline::csnd(real_t gm1) const
{
  return (gm1 >= gm1_low) ? csnd_rho(rho_gm1(gm1)) : poly.csnd(gm1);
}

real_t eos_barotr_spline::csnd_from_rho_gm1(real_t rho, real_t gm1) const 
{
  return (rho >= rho_low) ? csnd_rho(rho) : poly.csnd(gm1);
}

real_t eos_barotr_spline::temp(real_t gm1) const
{
  if (zerotemp) return 0.;
  return (gm1 <= gm1_low) ? temp0 : (*temp_gm1)(gm1);
}

real_t eos_barotr_spline::ye(real_t gm1) const
{
  if (!efrac_gm1)
    throw runtime_error("eos_barotr_table: electron fraction "
                        "not available.");
  return (gm1 >= gm1_low) ? (*efrac_gm1)(gm1) : efrac0;
}


auto eos_barotr_spline::descr_str() const -> std::string
{
  auto u = units_to_SI();
  std::ostringstream s;
  s.precision(15);
  s.setf(std::ios::scientific);
  s << "Interpolating spline EOS, "
    << "max. valid density =" 
    << (range_rho().max() * u.density())
    << " kg/m^3, "
    << "max. valid g-1 =" << range_gm1().max() << ", ";
    
    if (has_temp()) {
      if (is_zero_temp()) {
        s << "zero temperature";
      } else {
        s << "available";
      }
    } else {
      s << "not available";
    }
    s << ", ";
    if (is_isentropic()) {
      s << "isentropic, ";
    }
    
  s << "electron fraction " 
    << (has_efrac() ? "" : "not") << " available"
    << ". Below density of " 
    << (rho_low * u.density())
    << " kg/m^3"
    << " using " << poly.descr_str();
  return s.str();
}
  
auto EOS_Toolkit::make_eos_barotr_spline(
  func_t gm1_rho, func_t rho_gm1, func_t eps_gm1, func_t press_gm1, 
  func_t csnd_rho, func_t temp_gm1, func_t efrac_gm1, 
  bool isentropic, interval<real_t> rg_rho, real_t n_poly,
  units u, std::size_t pts_per_mag)
-> eos_barotr
{  
  using detail::interpol_logspl_impl;
  using detail::interpol_llogspl_impl;
  
  const size_t fac_pts_rho{ 5 };
  
  auto hm1_gm1{ 
    [&] (real_t gm1) -> real_t {
      return eps_gm1(gm1) + press_gm1(gm1) / rho_gm1(gm1);
    }
  };
  

  real_t rho_join{ rg_rho.min() };
  real_t gm1_join{ gm1_rho(rho_join) };
  real_t eps_join{ eps_gm1(gm1_join) };
  real_t p_join{ press_gm1(gm1_join) };
  real_t rhomax_poly{ 1.000001 * rho_join };
  auto poly = eos_barotr_gpoly::from_boundary(rho_join, eps_join, 
                               p_join, n_poly, rhomax_poly, u);

  const real_t gcorr{
    (poly.gm1_from_rho(rho_join) - gm1_join) / (1.0 + gm1_join)
  };
  auto gm1_new = [&] (real_t gm1o) -> real_t {
    return gm1o + gcorr * (1.0 + gm1o);
  };
  auto gm1_old = [&] (real_t gm1n) -> real_t {
    return gm1n - (gcorr / (1.0 + gcorr)) * (1.0 + gm1n);
  };

  const interval<real_t> rg_gm1{ gm1_new(gm1_join), 
                                 gm1_new(gm1_rho(rg_rho.max())) };

  if (rg_gm1.min()  <= 0) {
    throw std::range_error("eos_barotr_spline: invalid interval "
                           "requested for interpolation range");
  }
  
  auto get_npts = [&] (interval<real_t> rg) -> std::size_t {
    real_t lgr{ std::max(1.0, log10(rg.max() / rg.min())) };
    return pts_per_mag * std::size_t(lgr);
  };
  
  std::size_t npts_gm1{ get_npts(rg_gm1) };
  std::size_t npts_rho{ fac_pts_rho * npts_gm1 };
    
  auto sgm1{ 
    interpol_llogspl_impl::from_function(
      [&] (real_t rho) -> real_t {return gm1_new(gm1_rho(rho));},
      rg_rho, npts_rho)
  };
  
  auto srho{ 
    interpol_llogspl_impl::from_function(
      [&] (real_t gm1) -> real_t {return rho_gm1(gm1_old(gm1));},
      rg_gm1, npts_gm1)
  };
  
  auto seps{ 
    interpol_logspl_impl::from_function(
      [&] (real_t gm1) -> real_t {return eps_gm1(gm1_old(gm1));},
      rg_gm1, npts_gm1)
  };
  
  auto shm1{ 
    interpol_logspl_impl::from_function(
      [&] (real_t gm1) -> real_t {return hm1_gm1(gm1_old(gm1));},
      rg_gm1, npts_gm1)  
  };
  
  auto spress{ 
    interpol_llogspl_impl::from_function(
      [&] (real_t gm1) -> real_t {return press_gm1(gm1_old(gm1));},  
      rg_gm1, npts_gm1)
  };
  
  auto scsnd{ 
    interpol_logspl_impl::from_function(csnd_rho, rg_rho, npts_rho)
  };
  
  boost::optional<interpol_logspl_impl> stemp;
  if (temp_gm1) 
  {
    stemp = interpol_logspl_impl::from_function(
              [&] (real_t gm1) -> real_t {
                 return temp_gm1(gm1_old(gm1));
              },
              rg_gm1, npts_gm1);
  }
  
  boost::optional<interpol_logspl_impl> sefrac;
  if (efrac_gm1)  
  {
    sefrac = interpol_logspl_impl::from_function(
                 [&] (real_t gm1) {
                   return efrac_gm1(gm1_old(gm1));
                 },
                 rg_gm1, npts_gm1);
  }
  
  return eos_barotr{ 
    std::make_shared<eos_barotr_spline>(sgm1, srho, seps, 
       spress, shm1, scsnd, stemp, sefrac, isentropic, poly) 
  };
}


auto EOS_Toolkit::make_eos_barotr_spline(const eos_barotr& eos, 
              interval<real_t> rg_rho, real_t n_poly,
              std::size_t pts_per_mag)
-> eos_barotr
{
  func_t temp_gm1{nullptr};
  if (eos.has_temp())
  {
    temp_gm1 = [&eos](real_t gm1) {return eos.at_gm1(gm1).temp();};
  }
  func_t efrac_gm1{nullptr};
  if (eos.has_efrac())
  {
    efrac_gm1 = [&eos](real_t gm1) {return eos.at_gm1(gm1).ye();};
  }
  
  return make_eos_barotr_spline(
            [&eos](real_t rho) {return eos.at_rho(rho).gm1();},  
            [&eos](real_t gm1) {return eos.at_gm1(gm1).rho();},  
            [&eos](real_t gm1) {return eos.at_gm1(gm1).eps();},  
            [&eos](real_t gm1) {return eos.at_gm1(gm1).press();},  
            [&eos](real_t rho) {return eos.at_rho(rho).csnd();},  
            temp_gm1, efrac_gm1, 
            eos.is_isentropic(), rg_rho, n_poly,
            eos.units_to_SI(), pts_per_mag
         );
}

template<class T>
auto integrate_trapz(const std::vector<T>& x, const std::vector<T>& y, 
                     const T int_const=0)
-> std::vector<T>
{
  assert(x.size() == y.size());
  assert(x.size() > 0);
  
  std::vector<T> r(x.size());
  r[0] = int_const;
  for (std::size_t i=0; i+1<x.size(); ++i) 
  {
    r[i+1] = r[i] + (x[i+1] - x[i]) * (y[i] + y[i+1]) / 2;
  }
  return r;
}

auto compute_gm1_for_samples(const std::vector<real_t>& rho,
                             func_t eps_rho,
                             func_t press_rho)
-> std::vector<real_t>
{
  std::vector<real_t> vp, vy;
  for (auto r_ : rho) 
  {
    assert(r_ >0);
    const real_t p_{ press_rho(r_) };
    const real_t y_{ 1.0 / (r_ * (1.0 + eps_rho(r_)) + p_) };
    vp.push_back(p_);
    vy.push_back(y_);
  }
  
  auto gm1{ integrate_trapz(vp, vy) };
  
  for (real_t& x : gm1) x = expm1(x);
  
  return gm1;
}


auto ensure_resolution_log(const std::vector<real_t>& x, 
                           std::size_t ppmg) 
-> std::vector<real_t>
{
  std::vector<real_t> xf;

  const real_t mindlgx{ log(10.) / ppmg };
  for (std::size_t i=0; i+1<x.size(); ++i) 
  {    
    assert(x[i] > 0);
    assert(x[i+1] > x[i]);
    const real_t lgx0{ log(x[i]) };
    const real_t lgx1{ log(x[i+1]) };
    const real_t dlgx{ lgx1 - lgx0 };
    const int    nsub{ (int) ceil(dlgx / mindlgx) };
    const real_t sdlgx{ dlgx / nsub };
    for (int k=0; k<nsub; ++k)
    {
      xf.push_back(exp(lgx0 + k*sdlgx));
    }
  
  }
  xf.push_back(x.back());
  
  return xf;
}


auto EOS_Toolkit::make_eos_barotr_spline(
  const std::vector<real_t>& rho,
  const std::vector<real_t>& eps, const std::vector<real_t>& press, 
  const std::vector<real_t>& csnd, const std::vector<real_t>& temp, 
  const std::vector<real_t>& efrac, bool isentropic, 
  interval<real_t> rg_rho, real_t n_poly,
  units uc, std::size_t pts_per_mag)
-> eos_barotr
{  
  if (rho[0] <= 0) {
    throw std::runtime_error("Density must be strictly positive for"
    "sample points when constructing eos_barotr_spline EOS");
  }
  
  auto eps_rho{ make_interpol_pchip_spline(rho, eps) };
  auto press_rho{ make_interpol_pchip_spline(rho, press) };
  auto csnd_rho{ make_interpol_pchip_spline(rho, csnd) };

  auto rho_f{ ensure_resolution_log(rho, pts_per_mag) };
  
  auto gm1_f{ compute_gm1_for_samples(rho_f, eps_rho, press_rho) };

  auto rho_gm1{ make_interpol_pchip_spline(gm1_f, rho_f) };
  auto gm1_rho{ make_interpol_pchip_spline(rho_f, gm1_f) };

  auto eps_gm1 = [&] (real_t gm1) {return eps_rho(rho_gm1(gm1));};
  auto press_gm1 = [&] (real_t gm1) {return press_rho(rho_gm1(gm1));};

  auto temp_gm1 = [&]() -> func_t {
    if (temp.empty()) return nullptr;
    auto temp_rho{ make_interpol_pchip_spline(rho, temp) };
    return [=] (real_t gm1) {return temp_rho(rho_gm1(gm1));};
  }();

  auto efrac_gm1 = [&]() -> func_t {
    if (efrac.empty()) return nullptr;
    auto efrac_rho{ make_interpol_pchip_spline(rho, efrac) };
    return [=] (real_t gm1) {return efrac_rho(rho_gm1(gm1));};
  }();
  
  if (!gm1_rho.contains(rg_rho))
  {
    throw std::range_error("eos_barotr_spline: target density range "
                           "outside provided sample points");
  }
  
  return make_eos_barotr_spline(gm1_rho, rho_gm1, eps_gm1, 
                      press_gm1, csnd_rho, temp_gm1, efrac_gm1, 
                      isentropic, rg_rho, n_poly, uc, pts_per_mag);
}


auto compute_eps_for_samples(const std::vector<real_t>& rho,
                             func_t press_rho)
-> std::vector<real_t>
{
  std::vector<real_t> vy;
  for (auto r_ : rho) 
  {
    assert(r_ > 0);
    const real_t p_{ press_rho(r_) };
    vy.push_back( p_ / (r_ * r_) );
  }
  
  return integrate_trapz(rho, vy);
}


auto EOS_Toolkit::make_eos_barotr_spline(
  const std::vector<real_t>& rho,
  const std::vector<real_t>& press,
  const std::vector<real_t>& csnd,
  const std::vector<real_t>& temp,
  const std::vector<real_t>& efrac, 
  interval<real_t> rg_rho,
  real_t n_poly, real_t eps0, units uc,
  std::size_t pts_per_mag) 
-> eos_barotr
{
  if (rho[0] <= 0) {
    throw std::runtime_error("Density must be strictly positive for"
    "sample points when constructing eos_barotr_spline EOS");
  }
  
  auto press_rho{ make_interpol_pchip_spline(rho, press) };
  auto csnd_rho{ make_interpol_pchip_spline(rho, csnd) };

  auto rho_f{ ensure_resolution_log(rho, pts_per_mag) };
  
  
  const real_t eps_m{ n_poly * press_rho(rg_rho.min()) / rg_rho.min() };
  
  auto deps_f{ compute_eps_for_samples(rho_f, press_rho) };
  
  auto deps_rho{ make_interpol_pchip_spline(rho_f, deps_f) };
  auto eps_rho = [&] (real_t rho) {
    return deps_rho(rho) - deps_rho(rg_rho.min()) + eps_m;
  };
  
  auto gm1_f{ compute_gm1_for_samples(rho_f, eps_rho, press_rho) };

  auto rho_gm1{ make_interpol_pchip_spline(gm1_f, rho_f) };
  auto gm1_rho{ make_interpol_pchip_spline(rho_f, gm1_f) };

  auto eps_gm1 = [&] (real_t gm1) {return eps_rho(rho_gm1(gm1));};
  auto press_gm1 = [&] (real_t gm1) {return press_rho(rho_gm1(gm1));};

  auto temp_gm1 = [&]() -> func_t {
    if (temp.empty()) return nullptr;
    auto temp_rho{ make_interpol_pchip_spline(rho, temp) };
    return [=] (real_t gm1) {return temp_rho(rho_gm1(gm1));};
  }();

  auto efrac_gm1 = [&]() -> func_t {
    if (efrac.empty()) return nullptr;
    auto efrac_rho{ make_interpol_pchip_spline(rho, efrac) };
    return [=] (real_t gm1) {return efrac_rho(rho_gm1(gm1));};
  }();
  
  if (!gm1_rho.contains(rg_rho))
  {
    throw std::range_error("eos_barotr_spline: target density range "
                           "outside provided sample points");
  }
  
  return make_eos_barotr_spline(gm1_rho, rho_gm1, eps_gm1, 
                      press_gm1, csnd_rho, temp_gm1, efrac_gm1, 
                      true, rg_rho, n_poly, uc, pts_per_mag);
  
}

auto EOS_Toolkit::make_eos_barotr_spline(
  const std::vector<real_t>& gm1, const std::vector<real_t>& rho,
  const std::vector<real_t>& eps, const std::vector<real_t>& press, 
  const std::vector<real_t>& csnd, const std::vector<real_t>& temp, 
  const std::vector<real_t>& efrac, bool isentropic, 
  interval<real_t> rg_rho, real_t n_poly,
  units uc, std::size_t pts_per_mag)
-> eos_barotr
{  
  
  auto gm1_rho{ make_interpol_pchip_spline(rho, gm1) };
  auto rho_gm1{ make_interpol_pchip_spline(gm1, rho) };
  auto eps_gm1{ make_interpol_pchip_spline(gm1, eps) };
  auto press_gm1{ make_interpol_pchip_spline(gm1, press) };
  auto csnd_rho{ make_interpol_pchip_spline(rho, csnd) };
  
  func_t temp_gm1{nullptr};
  if (!temp.empty()) 
  {
    temp_gm1 = make_interpol_pchip_spline(gm1, temp);
  }

  func_t efrac_gm1{nullptr};
  if (!efrac.empty())
  {
    efrac_gm1 = make_interpol_pchip_spline(gm1, efrac);
  }
  
  if (!gm1_rho.contains(rg_rho))
  {
    throw std::range_error("eos_barotr_spline: target density range "
                           "outside provided sample points");
  }
  
  return make_eos_barotr_spline(gm1_rho, rho_gm1, eps_gm1, 
                      press_gm1, csnd_rho, temp_gm1, efrac_gm1, 
                      isentropic, rg_rho, n_poly, uc, pts_per_mag);
  
}

