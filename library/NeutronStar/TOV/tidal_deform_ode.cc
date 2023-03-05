#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include "tidal_deform_ode.h"
#include "solve_ode.h"

using namespace EOS_Toolkit;

namespace {

const EOS_Toolkit::real_t PI {
   boost::math::constants::pi<EOS_Toolkit::real_t>()
};

}



auto EOS_Toolkit::find_deform(eos_barotr eos, real_t gm1_center, 
            const std::vector<real_t>& dnu, 
            const std::vector<real_t>& rsqr, 
            const std::vector<real_t>& lambda, 
            real_t acc) -> spherical_star_tidal
{
    real_t gm1_switch{ gm1_center / 1.1 };
    real_t rho_switch{ eos.at_gm1(gm1_switch).rho() };
    tidal_ode tode(eos, gm1_center, dnu, rsqr, lambda, rho_switch);

    auto rtid{ integrate_ode_adapt(tode, acc, acc) };

    real_t dnu_switch{
      -std::log1p( (gm1_switch-gm1_center) / (1.+gm1_center) )
    };
    
    real_t y_switch{ rtid[tidal_ode::YM2] + 2. };
    
    tidal_ode2 tode2( eos, gm1_center, dnu, rsqr, lambda, 
                      dnu_switch, y_switch);

    auto rtid2{ integrate_ode_adapt(tode2, acc, acc) };

    return tode2.deformability(rtid2);
}



auto tidal_ode::gm1_from_dnu(real_t dnu) const -> real_t
{
  real_t gm1_raw{ gm1_center + (1.0 + gm1_center) * std::expm1(-dnu) };
  return std::max(gm1_raw, 0.0);
}

tidal_ode::tidal_ode(eos_barotr eos_, real_t gm1_center_, 
            const std::vector<real_t>& dnu_, 
            const std::vector<real_t>& rsqr_, 
            const std::vector<real_t>& lambda_,
            real_t rho_stop_)
: eos{eos_}, gm1_center{gm1_center_}, rho_stop{rho_stop_}
{
  if (!eos.is_isentropic()) {
    throw(std::runtime_error("Tidal deformability can only be"
                             "computed for isentropic EOS"));
  }
  
  std::vector<real_t> revrho, revlambda, revdnu, revrsqr, revmbr3;
  
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
    auto s{ eos.at_gm1(eos.range_gm1().limit_to(gm1)) };
    
    assert(s);
    assert(lambda >= 0.0);
    assert(rsqr >= 0.0);
    
    real_t rho{ s.rho() };
    real_t rho_e{ rho * (1.0 + s.eps()) };
    real_t mbr3{ m_by_r3(rsqr, lambda, rho_e) };
    
    revrho.push_back(rho);
    revlambda.push_back(lambda);
    revdnu.push_back(dnu);
    revrsqr.push_back(rsqr);
    revmbr3.push_back(mbr3);
  }
  
  
  dnu_rho = make_interpol_pchip_spline(revrho, revdnu);
  
  lambda_rho = make_interpol_pchip_spline(revrho, revlambda);

  rsqr_rho = make_interpol_pchip_spline(revrho, revrsqr);

  mbr3_rho = make_interpol_pchip_spline(revrho, revmbr3);


  assert(x_start()>x_end());
}



auto tidal_ode::m_by_r3(real_t rsqr, real_t lambda, 
                        real_t rho_e) const -> real_t
{
  if (rsqr <= 0) {
    return (4.0/3.0) * PI * rho_e;
  }
  return -0.5 * std::expm1(-2.0 * lambda) / rsqr;   
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
  assert(lambda_rho.range_x().contains(rho));
  real_t lambda{ std::max(0.0, lambda_rho(rho)) };
  real_t e2l{ std::exp(2.0 * lambda) };
  real_t rsqr{ rsqr_rho(rho) };
  real_t wtfac{ cs2 / rho };
  real_t mbyr3{ mbr3_rho(rho) };
  
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
  return std::max(gm1_raw, 0.0);
}



tidal_ode2::tidal_ode2(eos_barotr eos_, real_t gm1_center_, 
            const std::vector<real_t>& dnu_, 
            const std::vector<real_t>& rsqr_, 
            const std::vector<real_t>& lambda_, 
            real_t dnu0_, real_t y0_)
: eos{eos_}, gm1_center{gm1_center_}, dnu0{dnu0_}
{
  std::vector<real_t> rddy, rrho;
  assert(dnu_.size() == rsqr_.size());
  assert(dnu_.size() == lambda_.size());
  
  for (std::size_t k = dnu_.size()-1; k > 0; --k) 
  {
    real_t gm1{ gm1_from_dnu(dnu_[k]) };
    auto s{ eos.at_gm1(eos.range_gm1().limit_to(gm1)) };
    assert(s);
    
    rrho.push_back( s.rho() );
    
    real_t h{ 1. + s.hm1() };
    real_t p{ s.press() };
    real_t mbr3{ m_by_r3(rsqr_[k], lambda_[k]) };
  
    rddy.push_back( h / (p + mbr3 / (4*PI)) );
  }
  
  std::vector<real_t> rdy(rrho.size());
  rdy[0] = 0;
  for (size_t k=1; k<rrho.size(); ++k) {
    real_t drho{ (rrho[k]-rrho[k-1]) };
    assert(drho > 0);
    rdy[k] = rdy[k-1] + 0.5 * (rddy[k] + rddy[k-1]) * drho;
  }
  
  deltay_rho = make_interpol_pchip_spline(rrho, rdy);

  rsqr_dnu   = make_interpol_pchip_spline(dnu_, rsqr_);
  lambda_dnu = make_interpol_pchip_spline(dnu_, lambda_);

  real_t rho0{ 
    eos.at_gm1(eos.range_gm1().limit_to(gm1_from_dnu(dnu0))).rho() 
  }; 
  yhat0 = y0_ - deltay_rho(rho0);
}


auto tidal_ode2::dlnh_yhat(real_t dnu, real_t yhat) const -> real_t
{
  real_t gm1{ gm1_from_dnu(dnu) };
  real_t lambda{ lambda_dnu(dnu) };
  real_t rsqr{ rsqr_dnu(dnu) };
  
  //limit to range because to prevent roundoff errors causing trouble
  //when central density is at maximum of validity range
  auto s{ eos.at_gm1(eos.range_gm1().limit_to(gm1)) };
  assert(s);
  
  real_t rho{ s.rho() };
  real_t p{ s.press() };
  real_t eps{ s.eps() };
  real_t rho_e{ rho * (1. + eps) };
  real_t mbr3{ m_by_r3(rsqr, lambda) };
 
  real_t y{ yhat + deltay_rho(rho) };
  
  real_t a{ 4.0*PI * p + mbr3 };
  real_t b{ rsqr * std::exp(2.*lambda) };
  real_t c{ 2.*PI * (3.*rho_e + 11.*p) - 4.*mbr3 };
  real_t d{ (y + 3.) / (2. * b)  + mbr3 + 2.*PI * (p - rho_e) };
  
  return ((y-2.) * d + c) * 2. / a - 4. * b * a;  
}

void tidal_ode2::operator()(const state_t &s, state_t &dsdx, 
                  const real_t x) const
{
   dsdx[YHAT] = -dlnh_yhat(x, s[YHAT]);
}


auto tidal_ode2::initial_data() const -> state_t
{
  return {yhat0};  
}



auto tidal_ode2::deform_from_y_mbr(real_t y, real_t b) 
const -> spherical_star_tidal
{
  const double c0{ 
    2.0*b*(6.0 - 3.0*y + 3.0*b*(5.0*y-8.0) 
           + 2 * (b*b) * (13.0 - 11.0*y + b * (3.0*y-2.0) 
                          + 2.0 * (b*b) * (1.0+y)))
    +3.0 * pow((1.0 - 2.0*b), 2) 
         * (2.0 - y + 2.0*b*(y-1.0)) * log(1.0-2.0*b)
  };
  
  const double c1{ 
    ((8.0/5.0) * pow(b, 5)) * pow((1.0-2.0*b), 2) 
    * (2.0 - y + 2.0 * b * (y - 1.0))
  };
    
  real_t k2{ c1 / c0 };  
  real_t lt{ (2.0/3.0) * k2 / pow(b, 5) };
  
  return {k2, lt};
}


auto tidal_ode2::deformability(const state_t& surf) 
const -> spherical_star_tidal
{ 
  real_t dnu{ x_end() };
  real_t gm1{ gm1_from_dnu(dnu) };
  real_t lambda{ lambda_dnu(dnu) };
  
  auto s{ eos.at_gm1(eos.range_gm1().limit_to(gm1)) };
  assert(s);  
  real_t rho{ s.rho() };

  real_t mbr{ -0.5 * std::expm1(-2.0 * lambda) };
  real_t y{ surf[YHAT] + deltay_rho(rho) };
  
  return deform_from_y_mbr(y, mbr);
}




