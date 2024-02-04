#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include "tidal_deform_ode.h"
#include "solve_ode.h"
#include "integration.h"

using namespace EOS_Toolkit;

namespace {

const EOS_Toolkit::real_t PI {
   boost::math::constants::pi<EOS_Toolkit::real_t>()
};

}


auto EOS_Toolkit::k2_from_ym2_mbr_stable(real_t ym2, real_t mbr, 
                              real_t b_thresh)  -> real_t
{
  return  tidal_ode2::k2_from_ym2_mbr_interp(ym2, mbr, b_thresh); 
}

auto tidal_ode::gm1_from_dnu(real_t dnu) const -> real_t
{
  real_t gm1_raw{ gm1_center + (1.0 + gm1_center) * std::expm1(-dnu) };
  return eos.range_gm1().limit_to(gm1_raw);
}

tidal_ode::tidal_ode(eos_barotr eos_, 
            const spherical_star_info& prop_,
            const std::vector<real_t>& dnu_, 
            const std::vector<real_t>& rsqr_, 
            const std::vector<real_t>& lambda_,
            real_t rho_stop_)
: eos{eos_}, gm1_center{prop_.center_gm1}, 
  rho_start{prop_.center_rho}, rho_stop{rho_stop_}

{
  if (!eos.is_isentropic()) {
    throw(std::runtime_error("Tidal deformability can only be"
                             "computed for isentropic EOS"));
  }
  assert(rsqr_[1]>0);
  std::vector<real_t> revgm1, revlambda, revrsqr, revmbr3;
  
  auto ilambda = lambda_.rbegin();
  auto irsqr = rsqr_.rbegin();
  for (auto idnu = dnu_.rbegin(); idnu != dnu_.rend(); ++idnu) 
  {
    assert(ilambda != lambda_.rend());
    assert(irsqr != rsqr_.rend());
    
    real_t dnu{ *idnu };
    real_t lambda{ *ilambda++ };
    real_t rsqr{ *irsqr++ };

    real_t gm1{ gm1_from_dnu(dnu) };
    auto s{ eos.at_gm1(gm1) };
    
    assert(s);
    assert(lambda >= 0.0);
    assert(rsqr >= 0.0);
    
    real_t rho{ s.rho() };
    real_t rho_e{ rho * (1.0 + s.eps()) };
    real_t mbr3{ 
      rsqr >= rsqr_[1] ? m_by_r3(rsqr, lambda, rho_e) 
                       : m_by_r3_origin(rho_e) 
    };
    
    revgm1.push_back(gm1);
    revlambda.push_back(lambda);
    revrsqr.push_back(rsqr);
    revmbr3.push_back(mbr3);
  }
  
  
  lambda_gm1 = make_interpol_pchip_spline(revgm1, revlambda);

  rsqr_gm1 = make_interpol_pchip_spline(revgm1, revrsqr);

  mbr3_gm1 = make_interpol_pchip_spline(revgm1, revmbr3);


  assert(x_start()>x_end());
}



auto tidal_ode::m_by_r3(real_t rsqr, real_t lambda, 
                        real_t rho_e) const -> real_t
{
  assert(rsqr > 0);
  return -0.5 * std::expm1(-2.0 * lambda) / rsqr;   
}


auto tidal_ode::m_by_r3_origin(real_t rho_e) -> real_t
{
  return (4.0/3.0) * PI * rho_e;
}

auto tidal_ode::drho_y(real_t rho_, real_t ym2) const -> real_t
{
  assert(rho_>0);
  const real_t rho{ eos.range_rho().limit_to(rho_) };
  
  auto s{ eos.at_rho(rho) };
  assert(s);
  real_t h{ s.hm1() + 1. }; 
  real_t p{ s.press() }; 
  real_t eps{ s.eps() }; 
  real_t cs2{ std::pow(s.csnd(), 2) }; 
   
  real_t rho_e{ rho * (1.0 + eps) };
  real_t gm1{ lambda_gm1.range_x().limit_to(s.gm1()) }; 
  real_t lambda{ std::max(0.0, lambda_gm1(gm1)) };
  real_t e2l{ std::exp(2.0 * lambda) };
  real_t rsqr{ rsqr_gm1(gm1) };
  real_t wtfac{ cs2 / rho };
  real_t mbyr3{ mbr3_gm1(gm1) };
  
  real_t a{ 4.0*PI * p + mbyr3 };
  real_t b{ 2.0 * rsqr * (mbyr3 + 2.0*PI * (p - rho_e)) };
  real_t c{ 4.0*PI * (3.0 * rho_e + 11.0 * p) - 8.0 * mbyr3 };
  real_t d{ 4.0 * rsqr * wtfac * e2l * a };
  real_t g{ (ym2+2.0 + 3.0) / e2l + b };
  
  
  real_t f{0.0};
  if (rsqr>0) {
    f = wtfac * ym2 / rsqr;
  }
  else {    
    f = (-4.0*PI/7.0) * (
                (11.0 * h - (32.0/3.0) * (1.0 + eps) ) * cs2 + h);
  }


  real_t res{ (4.0*PI * h +  (f * g + wtfac * c)) / a - d };
  assert(std::isfinite(res));
  return res;
}

auto tidal_ode::initial_data() const -> state_t
{
  return {0.0};
}


void tidal_ode::operator()(const state_t &s , state_t &dsdx, 
                           const real_t x) const
{
  assert(std::isfinite(s[YM2]));
  dsdx[YM2] = drho_y(x, s[YM2]);
  assert(std::isfinite(dsdx[YM2]));
}



auto tidal_ode2::m_by_r3(real_t rsqr, real_t lambda) const -> real_t
{
  assert(rsqr>0);
  return -0.5 * std::expm1(-2.0 * lambda) / rsqr;   
}


auto tidal_ode2::gm1_from_dnu(real_t dnu) const -> real_t
{
  real_t gm1_raw{ gm1_center + (1.0 + gm1_center) * std::expm1(-dnu) };
  //limit to range because to prevent roundoff errors causing trouble
  //when central density is at maximum of validity range and at the 
  //surface
  return eos.range_gm1().limit_to(gm1_raw);
}




tidal_ode2::tidal_ode2(eos_barotr eos_, 
            const spherical_star_info& prop,
            const std::vector<real_t>& dnu_, 
            const std::vector<real_t>& rsqr_, 
            const std::vector<real_t>& lambda_, 
            real_t rho0_, real_t z0_)
: eos{eos_}, gm1_center{prop.center_gm1}
{
  assert(dnu_.size() == rsqr_.size());
  assert(dnu_.size() == lambda_.size());

  if (!eos.is_isentropic()) {
    throw(std::runtime_error("Tidal deformability can only be"
                             "computed for isentropic EOS"));
  }

  const real_t gm10{ eos.at_rho(rho0_).gm1() };
  dnu0 = -std::log1p((gm10 - gm1_center) / (1. + gm1_center));

  rsqr_dnu   = make_interpol_pchip_spline(dnu_, rsqr_);
  lambda_dnu = make_interpol_pchip_spline(dnu_, lambda_);

  std::vector<real_t> t_gm1, t_mbr3;
  for (std::size_t k = dnu_.size()-1; k > 0; --k) 
  {
    const real_t gm1{ gm1_from_dnu(dnu_[k]) };
    const real_t mbr3{ m_by_r3(rsqr_[k], lambda_[k]) };
    
    t_gm1.push_back( gm1 );
    t_mbr3.push_back(mbr3);
  }
  auto mbr3_gm1 { make_interpol_pchip_spline(t_gm1, t_mbr3) };

  std::vector<real_t> rddy, rrho;

  real_t dlgrho_max { 10. / dnu_.size() };
  
  rrho.push_back(eos.at_gm1(t_gm1[0]).rho());
  for (std::size_t j=1; j<dnu_.size()-2; ++j) 
  {
    real_t rhoa { eos.at_gm1(t_gm1[j]).rho() };
    rrho.push_back( rhoa );
    
    real_t rhob { eos.at_gm1(t_gm1[j+1]).rho() };
    real_t lgrhoa{ log(rhoa) };
    real_t lgrhob{ log(rhob) };
    real_t dlgrho { lgrhob - lgrhoa };
    
    
    if (dlgrho > dlgrho_max)
    {
      const real_t n{ ceil(dlgrho/dlgrho_max) };
      for (int i=1; i<n; ++i) 
      {
        rrho.push_back( exp(lgrhoa + dlgrho * i / real_t(n)) );
      }
    }
  }
  rrho.push_back(eos.at_gm1(t_gm1[dnu_.size()-2]).rho()  );
  
  for (auto rho :  rrho) 
  {
      auto s{ eos.at_rho(rho) };
      assert(s);
      real_t h{ 1. + s.hm1() };
      real_t p{ s.press() };
      real_t gm1{ mbr3_gm1.range_x().limit_to(s.gm1()) };
      
      rddy.push_back( h / (p + mbr3_gm1(gm1) / (4.*PI)) );
  }
  auto rdy = integrate_order3(rrho, rddy);    
  deltay_rho = make_interpol_pchip_spline(rrho, rdy);

   
  zhat0 = z0_ - deltay_rho(rho0_);
}


auto tidal_ode2::dlnh_zhat(real_t dnu, real_t zhat) const -> real_t
{
  real_t gm1{ gm1_from_dnu(dnu) };
  real_t lambda{ lambda_dnu(dnu) };
  real_t rsqr{ rsqr_dnu(dnu) };
  

  auto s{ eos.at_gm1(gm1) };
  assert(s);
  
  real_t rho{ s.rho() };
  real_t p{ s.press() };
  real_t eps{ s.eps() };
  real_t rho_e{ rho * (1. + eps) };
  real_t mbr3{ m_by_r3(rsqr, lambda) };
 
  real_t z{ zhat + deltay_rho(rho) };
  
  real_t a{ 4.0*PI * p + mbr3 };
  real_t b{ rsqr * std::exp(2.*lambda) };
  real_t c{ 2.*PI * (3.*rho_e + 11.*p) - 4.*mbr3 };
  real_t d{ (z + 5.) / (2. * b)  + mbr3 + 2.*PI * (p - rho_e) };
  
  return (z * d + c) * 2. / a - 4. * b * a;  
}

void tidal_ode2::operator()(const state_t &s, state_t &dsdx, 
                  const real_t x) const
{
   dsdx[ZHAT] = -dlnh_zhat(x, s[ZHAT]);
}


auto tidal_ode2::initial_data() const -> state_t
{
  return {zhat0};  
}



auto tidal_ode2::k2_from_ym2_mbr_taylor_exp(real_t z, real_t b) 
-> real_t
{
    const real_t zp5{ z + 5. };
    const real_t p0{ -z / 2. };
    const real_t p1{ 
      (z * (z + 10.) + 10.) / 2.
    };
    const real_t p2{ 
      (((z + 40.) * z + 140.)*z + 140.) / 14.
    };
    const real_t p3{ 
      ((((z  + 70. ) * z + 480.) * z + 1160.) * z + 980.) / 14.
    };
    const real_t p4{ 
      (((((5.*z + 540.) * z + 5592.) * z + 22876.) * z
       + 42756.) * z + 30912.) * 5. / 294.
    };
    const real_t p5{ 
      ((((((33. * z + 5090.) * z + 71020.) * z + 416280.) * z 
           + 1249720.) * z + 1918280.) * z + 1207920.) / 294.
    };
    const real_t c{ b/zp5 };
    const real_t p{ 
      (p0 + c * (p1 + c * (p2 + c * (p3 + c * (p4 + c * p5))))) / zp5
    };
    return  p * pow(1. - 2. * b, 2);
}

namespace {
auto lgsigma(real_t x) -> real_t 
{
  const real_t e{ expm1(x) };
  return log1p((x - e) / (1. + e));
}
}

  
auto tidal_ode2::k2_from_ym2_mbr(real_t z, real_t b) -> real_t
{
  assert(b > 0);
  const real_t s { lgsigma(-2.*b) / (-2.*b*b) };
  const real_t c1 {  
    -(4./5.) * pow(b, 3) * (2.*b*(z+1.) - z)
  };  
  const real_t c0 { 
    b * (b * (b * (24.*s*(z+1.) - 4.*z - 12.) 
              + ((-36.*z -24.) * s + 18.*z + 16.))
         + ((18.*z + 6.) * s -14.*z - 6.)) + 3.*z*(1.-s)
  };
  return  pow(1.-2.*b, 2) * c1 / c0;
}

auto tidal_ode2::k2_from_ym2_mbr_interp(real_t z, real_t b, 
                                     real_t b_thresh) -> real_t
{   
    const real_t f { k2_from_ym2_mbr(z, b) };
    
    if (b > b_thresh) return f;
    
    const real_t a { k2_from_ym2_mbr_taylor_exp(z, b) };
    
    const real_t e { 
       k2_from_ym2_mbr(z, b_thresh) 
       - k2_from_ym2_mbr_taylor_exp(z, b_thresh)
    };
    
    return a + e*pow(b/b_thresh,6);
}



auto tidal_ode2::deform_from_ym2_mbr(real_t ym2, real_t b) 
-> spherical_star_tidal
{
  const real_t k2{ k2_from_ym2_mbr_interp(ym2,b) };
    
  const real_t lt{ (2.0/3.0) * k2 / pow(b, 5) };
  
  return {k2, lt};
}


auto tidal_ode2::deformability(const state_t& surf) 
const -> spherical_star_tidal
{ 
  real_t dnu{ x_end() };
  real_t lambda{ lambda_dnu(dnu) };
  
  real_t mbr{ -0.5 * std::expm1(-2.0 * lambda) };
  real_t z{ surf[ZHAT] };
  
  return deform_from_ym2_mbr(z, mbr);
}




