#include <cassert>
#include <algorithm>
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include "tov_ode.h"


using namespace EOS_Toolkit;

namespace {

const EOS_Toolkit::real_t PI {
   boost::math::constants::pi<EOS_Toolkit::real_t>()
};


}
  
tov_ode::tov_ode(real_t rho_center_, eos_barotr eos_)
: eos{std::move(eos_)}
{
  auto e{ eos.at_rho(rho_center_) };
  if (!e) {
    throw std::runtime_error("TOV central density outside EOS range");
  }
  gm1_center = e.gm1();
  hm1_center = e.hm1();
  rho_center = e.rho();
  rsqr_norm  = std::log1p(gm1_center) / rho_center;
}


auto tov_ode::x_end() const -> real_t
{
  real_t x{ std::log1p(gm1_center) };
  assert(x>0);
  return x;
}

auto tov_ode::gm1_from_x(real_t x) const -> real_t
{
  real_t gm1_raw{ gm1_center + (1.0 + gm1_center) * std::expm1(-x) };
  return std::max(gm1_raw, 0.0);
}

auto tov_ode::m_by_r3(real_t rsqr, real_t lambda, 
                      real_t rho_e) -> real_t
{
  assert(rsqr >= 0);
  if (rsqr == 0) {
    return (4.0/3.0) * PI * rho_e;
  }
  return -0.5 * std::expm1(-2.0 * lambda) / rsqr;   
}

auto tov_ode::grav_mass(real_t r, real_t lambda) -> real_t
{
  return -0.5 * r * std::expm1(-2.0 * lambda);
}

auto tov_ode::bary_mass(real_t r, real_t lambda, real_t ybnd) -> real_t
{
  return grav_mass(r, lambda) + ybnd * r;
}

auto tov_ode::prop_vol(real_t r, real_t yvol) -> real_t
{
  return yvol * r;
}

auto tov_ode::moment_inertia(real_t r_s, real_t omega1_s, 
                             real_t omega2_s) -> real_t
{
  return std::pow(r_s, 3) / (3.0 * omega1_s / omega2_s + 2.0);
}

auto tov_ode::dx_rsqr(real_t lambda, real_t press,
                      real_t mbyr3) -> real_t
{
  return 2.0 / (std::exp(2.0 * lambda) * (4.0*PI * press + mbyr3));
}

auto tov_ode::dx_lambda(real_t press, real_t rho_e, 
                        real_t mbyr3) -> real_t
{
  return (4.0*PI * rho_e - mbyr3) / (4.0*PI * press + mbyr3);
}



auto tov_ode::ebnd_by_r3(real_t ybnd, real_t rsqr, real_t rho, 
                         real_t eps) -> real_t
{
  assert(rsqr >= 0);
  if (rsqr == 0) {
    return (-4.0/3.0) * PI * rho * eps;
  }
  return ybnd / rsqr;   
}

auto tov_ode::drsqr_ybnd(real_t ybnd, real_t rsqr, real_t lambda, 
                         real_t rho, real_t eps) -> real_t
{
  const real_t ebyr3{ ebnd_by_r3(ybnd, rsqr, rho, eps) }; 
  
  return 2*PI * rho * (std::expm1(lambda) - eps) - ebyr3 / 2.0;
}

auto tov_ode::vol_by_r3(real_t yvol, real_t rsqr) -> real_t
{
  assert(rsqr >= 0);
  if (rsqr == 0) {
    return (4.0/3.0) * PI;
  }
  return yvol / rsqr;   
}

auto tov_ode::drsqr_yvol(real_t yvol, real_t rsqr, 
                         real_t lambda) -> real_t
{
  const real_t vbyr3{ vol_by_r3(yvol, rsqr) }; 
  
  return 2.0*PI * std::exp(lambda)  - vbyr3 / 2.0;
}

auto tov_ode::drsqr_omega2(real_t rsqr, real_t lambda, real_t rho, 
                           real_t hm1, real_t drsqr_omega1, 
                           real_t omega1) -> real_t
{
  const real_t a{ 4.0*PI * rho * (1.+hm1) * std::exp(2.0*lambda) };
  return 0.5 * (rsqr * a - 3.0) * drsqr_omega1 + a * omega1;
}


auto tov_ode::drsqr_omega1( real_t rsqr, real_t omega2) const 
-> real_t
{
  assert(rsqr >= 0);
  if (rsqr == 0) {
    return (8./5.) * PI * rho_center * (1. + hm1_center) * 1.0;
  }
  return omega2 / rsqr;   
}


auto tov_ode::initial_data() const -> state_t
{
  state_t s;
  s[RSQR]   = 0.;
  s[LAMBDA] = 0.;
  s[YBND]   = 0.;
  s[YVOL]   = 0.;
  s[OMEGA1] = 1.0;
  s[OMEGA2] = 0.;
  return s;
}



void tov_ode::operator()(const state_t &s , state_t &dsdx, 
                         const real_t x) const
{
  //limit to range because to prevent roundoff errors causing trouble
  //when central density is at maximum of validity range
  auto e{ eos.at_gm1(eos.range_gm1().limit_to(gm1_from_x(x))) };
  assert(e);
  
  const real_t press{ e.press() };
  const real_t eps{ e.eps() };
  const real_t rho{ e.rho() };
  const real_t hm1{ e.hm1() };
  const real_t rho_e{ (eps  + 1.0) * rho };
  const real_t rsqr{ s[RSQR] * rsqr_norm };
  assert(s[RSQR] >= 0);
  assert(rsqr >= 0);
  const real_t mbyr3{ m_by_r3(rsqr, s[LAMBDA], rho_e) };
  const real_t volbyr{ s[YVOL] * rsqr_norm };
  const real_t dr2_w1{ drsqr_omega1(rsqr, s[OMEGA2] / rsqr_norm) }; 
  const real_t dx_r2{ dx_rsqr(s[LAMBDA], press, mbyr3) };
  dsdx[LAMBDA] = dx_lambda(press, rho_e, mbyr3);
  dsdx[RSQR]   = dx_r2 / rsqr_norm;
  assert(dsdx[RSQR] >= 0);
  dsdx[YBND]   = dx_r2 
                  * drsqr_ybnd(s[YBND], rsqr, s[LAMBDA], rho, eps); 
  dsdx[YVOL]   = dsdx[RSQR] * drsqr_yvol(volbyr, rsqr, s[LAMBDA]);
  dsdx[OMEGA1] = dx_r2 * dr2_w1;
  dsdx[OMEGA2] = rsqr_norm * dx_r2 
                  * drsqr_omega2(rsqr, s[LAMBDA], rho, 
                                 hm1, dr2_w1, s[OMEGA1]);
      
}


auto tov_ode::star(const state_t& surf) const 
-> spherical_star_info
{
  const real_t rc{ std::sqrt(surf[RSQR] * rsqr_norm) };
  const real_t surf_nu{ - surf[LAMBDA] };
  const real_t center_nu{ surf_nu - x_end() };
  const real_t mg{ grav_mass(rc, surf[LAMBDA]) };
  const real_t ebnd{ surf[YBND] * rc };
  const real_t pvol{ prop_vol(rc, surf[YVOL] * rsqr_norm) };
  const real_t mi{ moment_inertia(rc, surf[OMEGA1], 
                                  surf[OMEGA2] / rsqr_norm) };
  
  return {rho_center, gm1_center, center_nu, mg, ebnd, rc, pvol, mi};
}


void tov_ode::observer::operator()(const state_t& snew, value_t xnew) 
{
  dnu.push_back(xnew);
  rsqr.push_back(snew[tov_ode::RSQR]* rsqr_norm);
  lambda.push_back(snew[tov_ode::LAMBDA]);
  ebnd_by_r.push_back(snew[tov_ode::YBND]);
  pvol_by_r.push_back(snew[tov_ode::YVOL] * rsqr_norm);
}
