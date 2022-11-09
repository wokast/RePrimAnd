#include <stdexcept>
#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/roots.hpp>
#include "spherical_stars_internals.h"
#include "tov_ode.h"

namespace {

const EOS_Toolkit::real_t PI {
   boost::math::constants::pi<EOS_Toolkit::real_t>()
};

}

namespace EOS_Toolkit {


auto g_rcrc_from_lambda(real_t lambda) -> real_t
{
  return std::exp(2.0*lambda);  
}

auto g_phph_from_rc_th(real_t rc, real_t th) -> real_t
{
  return std::pow(rc * std::sin(th), 2);  
}

auto g_thth_from_rc(real_t rc) -> real_t
{
  return rc*rc;
}

auto lapse_from_nu(real_t nu) -> real_t
{
  return std::exp(nu);
}


spherical_star_profile::spherical_star_profile(eos_barotr eos_,
                                               real_t surf_radius_)
: _eos{std::move(eos_)}, _surf_radius{surf_radius_}
{}

auto spherical_star_profile::eos() const -> const eos_barotr&
{
  return _eos;
}


auto spherical_star_profile::surf_circ_radius() const -> real_t
{
  return _surf_radius;
}


auto spherical_star_profile::state_from_rc(real_t rc) const 
-> eos_barotr::state
{
  return eos().at_gm1(gm1_from_rc(rc));
}


namespace details {


tov_profile::tov_profile(eos_barotr eos_, const spherical_star_info &p_, 
                         vec_t rsqr_, vec_t delta_nu_, vec_t lambda_,
                         vec_t ybnd_, vec_t yvol_)
: spherical_star_profile{std::move(eos_), std::sqrt(rsqr_.back())}, 
  lambda_rsqr{make_interpol_pchip_spline(rsqr_, lambda_)}, 
  delta_nu_rsqr{make_interpol_pchip_spline(rsqr_, delta_nu_)}, 
  ybnd_rsqr{make_interpol_pchip_spline(rsqr_, ybnd_)}, 
  yvol_rsqr{make_interpol_pchip_spline(rsqr_, yvol_)}, 
  gm1_c{p_.center_gm1}, nu_c{p_.center_nu}, 
  mgrav{p_.grav_mass}, ebind{p_.binding_energy}
{}



void tov_profile::validate_rc(real_t rc) const
{
  if (rc<0) { 
    throw std::runtime_error(
                 "evaluating star profile at negative radius");
  }
}


auto tov_profile::center_gm1() const -> real_t 
{
  return gm1_c;
}

auto tov_profile::nu_from_rc_outside(real_t rc) const -> real_t
{
  return 0.5 * std::log1p( -2.0 * mgrav / rc );
}

auto tov_profile::nu_from_rc(real_t rc) const -> real_t
{
  validate_rc(rc);
  if (rc >= surf_circ_radius()) 
  {
    return nu_from_rc_outside(rc);  
  }
  return nu_c + delta_nu_rsqr(rc*rc);
}

auto tov_profile::lambda_from_rc(real_t rc) const -> real_t
{
  validate_rc(rc);
  if (rc >= surf_circ_radius()) 
  {
    return -nu_from_rc_outside(rc);  
  }
  return lambda_rsqr(rc*rc);
}

auto tov_profile::gm1_from_rc(real_t rc) const -> real_t
{
  validate_rc(rc);
  if (rc >= surf_circ_radius()) 
  {
    return 0.;  
  }
  const real_t delta_nu{ delta_nu_rsqr(rc*rc) };
  const real_t gm1_raw { 
    gm1_c + (1.0 + gm1_c) * std::expm1(-delta_nu) 
  };
  return std::max(gm1_raw, 0.0);
}
  
auto tov_profile::pvol_vacuum(real_t rc) const -> real_t
{
  const real_t rsqr{ rc*rc };  
  const real_t s{ rc * sqrt(1.0 - 2.0 * mgrav / rc) };
  const real_t a{ 15 * std::pow(mgrav, 3) * log(s + rc - mgrav) };
  const real_t b{ s * (mgrav * (15 * mgrav + 5 * rc) + 2 * rsqr) };
  return (a + b) * 4.0*PI / 6.0;  
}

  
auto tov_profile::pvol_from_rc(real_t rc) const -> real_t
{
  validate_rc(rc);
  const real_t rs{ surf_circ_radius() };
  if (rc > rs) 
  {
    return rs * yvol_rsqr(rs*rs) + pvol_vacuum(rc) - pvol_vacuum(rs);
  }
  return rc * yvol_rsqr(rc*rc);
}


auto tov_profile::mbary_from_rc(real_t rc) const -> real_t
{
  validate_rc(rc);
  if (rc >= surf_circ_radius()) 
  {
    return mgrav + ebind;
  }
  const real_t rsqr{ rc*rc };
  const real_t mg{ -0.5 * rc * std::expm1(-2.0 * lambda_rsqr(rsqr)) };
  return mg + rc * ybnd_rsqr(rsqr);
}
  

  
} // namespace details

} // namespace EOS_Toolkit

