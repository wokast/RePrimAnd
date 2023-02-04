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
  const real_t rho_stop;
  
  interpolator dnu_rho;
  interpolator lambda_rho;
  interpolator rsqr_rho;
  interpolator mbr3_rho;
  
  auto gm1_from_dnu(real_t dnu) const -> real_t;

  auto m_by_r3(real_t rsqr, real_t lambda, 
               real_t rho_e) const -> real_t;

  auto drho_y(real_t rho, real_t y) const -> real_t;

  public:
  
  enum index_t {YM2=0};

  using state_t = std::array<real_t, 1>;

  tidal_ode(eos_barotr eos_, real_t gm1_center_, 
            const std::vector<real_t>& dnu, 
            const std::vector<real_t>& rsqr, 
            const std::vector<real_t>& lambda,
            real_t rho_stop_);

  auto x_start() const -> real_t {return lambda_rho.range_x().max();}
  auto x_end() const -> real_t {return rho_stop;}
      
  void operator()(const state_t &s , state_t &dsdx, 
                  const real_t x) const;

  auto initial_data() const -> state_t;
  
};


class tidal_ode2 {
  const eos_barotr eos;
  
  const real_t gm1_center;
  const real_t dnu0;
  real_t yhat0;
  
  interpolator deltay_rho;
  interpolator rsqr_dnu;
  interpolator lambda_dnu;
  
  auto gm1_from_dnu(real_t dnu) const -> real_t;

  auto m_by_r3(real_t rsqr, real_t lambda) const -> real_t;

  auto dlnh_yhat(real_t dnu, real_t y) const -> real_t;

  auto deform_from_y_mbr(real_t y, real_t b) const -> spherical_star_tidal;

  public:
  
  enum index_t {YHAT=0};

  using state_t = std::array<real_t, 1>;

  tidal_ode2(eos_barotr eos_, real_t gm1_center_, 
            const std::vector<real_t>& dnu, 
            const std::vector<real_t>& rsqr, 
            const std::vector<real_t>& lambda, 
            real_t dnu0_, real_t y0_);

  auto x_start() const -> real_t {return dnu0;}
  auto x_end() const -> real_t {return rsqr_dnu.range_x().max();}
      
  void operator()(const state_t &s , state_t &dsdx, 
                  const real_t x) const;

  auto initial_data() const -> state_t;
  
  auto deformability(const state_t& surf) const -> spherical_star_tidal; 
};


auto find_deform(eos_barotr eos_, real_t gm1_center, 
            const std::vector<real_t>& dnu, 
            const std::vector<real_t>& rsqr, 
            const std::vector<real_t>& lambda, 
            real_t acc) -> spherical_star_tidal;

}

#endif
