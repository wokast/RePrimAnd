#ifndef TIDAL_DEFORM_ODE_H
#define TIDAL_DEFORM_ODE_H
#include <vector>
#include "config.h"
#include "interpol.h"
#include "tov_ode.h"
#include "spherical_stars.h"

namespace EOS_Toolkit {

class tidal_ode {
  const eos_barotr eos;
  
  const real_t gm1_center;
  const real_t rho_start;
  const real_t rho_stop;
  
  interpolator lambda_gm1;
  interpolator rsqr_gm1;
  interpolator mbr3_gm1;
  
  auto gm1_from_dnu(real_t dnu) const -> real_t;

  auto m_by_r3(real_t rsqr, real_t lambda, 
               real_t rho_e) const -> real_t;
  static auto m_by_r3_origin(real_t rho_e) -> real_t;

  auto drho_y(real_t rho, real_t y) const -> real_t;

  public:
  using value_t = real_t;
  enum index_t {YM2=0};

  using state_t = std::array<real_t, 1>;

  tidal_ode(eos_barotr eos_, 
            const spherical_star_info& prop,
            const std::vector<real_t>& dnu, 
            const std::vector<real_t>& rsqr, 
            const std::vector<real_t>& lambda,
            real_t rho_stop_);

  auto x_start() const -> real_t {return rho_start;}
  auto x_end() const -> real_t {return rho_stop;}
      
  void operator()(const state_t &s , state_t &dsdx, 
                  const real_t x) const;

  auto initial_data() const -> state_t;
  
};


class tidal_ode2 {
  const eos_barotr eos;
  
  const real_t gm1_center;
  real_t dnu0;
  real_t zhat0;
  
  interpolator deltay_rho;
  interpolator rsqr_dnu;
  interpolator lambda_dnu;
  
  auto gm1_from_dnu(real_t dnu) const -> real_t;

  auto m_by_r3(real_t rsqr, real_t lambda) const -> real_t;

  auto dlnh_zhat(real_t dnu, real_t zhat) const -> real_t;

  static auto k2_from_ym2_mbr(real_t ym2, real_t b) -> real_t;
  
  static auto k2_from_ym2_mbr_taylor_exp(real_t ym2, real_t b) 
  -> real_t;
  
  
  static auto deform_from_ym2_mbr(real_t ym2, real_t b) 
  -> spherical_star_tidal;

  public:
  using value_t = real_t;
  enum index_t {ZHAT=0};

  using state_t = std::array<real_t, 1>;

  tidal_ode2(eos_barotr eos_, 
            const spherical_star_info& prop,
            const std::vector<real_t>& dnu, 
            const std::vector<real_t>& rsqr, 
            const std::vector<real_t>& lambda, 
            real_t rho0_, real_t z0_);

  auto x_start() const -> real_t {return dnu0;}
  auto x_end() const -> real_t {return rsqr_dnu.range_x().max();}
      
  void operator()(const state_t &s , state_t &dsdx, 
                  const real_t x) const;

  auto initial_data() const -> state_t;
  
  auto deformability(const state_t& surf) const -> spherical_star_tidal; 

  static auto k2_from_ym2_mbr_interp(real_t ym2, real_t b, 
                              real_t b_thresh=5e-2) -> real_t;

};


}

#endif
