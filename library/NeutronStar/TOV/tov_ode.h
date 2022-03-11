#ifndef TOV_ODE_H
#define TOV_ODE_H
#include <vector>
#include "config.h"
#include "eos_barotropic.h"
#include "spherical_stars.h"

namespace EOS_Toolkit {


class tov_ode {
  eos_barotr eos;
  real_t gm1_center;
  real_t hm1_center;
  real_t rho_center;
  real_t rsqr_norm;
  
  public:
  enum index_t {RSQR=0, LAMBDA=1, YBND=2, YVOL=3, OMEGA1=4,OMEGA2=5};

  using value_t = real_t;
  using state_t = std::array<value_t, 6>;
  
  class observer {
    public:
    using state_t = tov_ode::state_t;
    using value_t = tov_ode::value_t;
    
    std::vector<value_t> dnu, rsqr, lambda, ebnd_by_r, pvol_by_r;
    
    observer(const tov_ode& ode) : rsqr_norm{ode.rsqr_norm} {}

    void operator()(const state_t& snew, value_t xnew);
    
    private:
    value_t rsqr_norm;
  };
  

  tov_ode(real_t rho_center_, eos_barotr eos_);

  auto center_gm1() const -> real_t {return gm1_center;}
  auto gm1_from_x(real_t x) const -> real_t;
  auto x_start() const -> real_t {return 0.0;}
  auto x_end() const -> real_t;
  
  static auto m_by_r3(real_t rsrq, real_t lambda, 
                      real_t rho_e) -> real_t;
  
  static auto grav_mass(real_t r, real_t lambda) -> real_t;
  static auto bary_mass(real_t r, real_t lambda, real_t ybnd) -> real_t;
  static auto prop_vol(real_t r, real_t yvol) -> real_t;
  static auto moment_inertia(real_t r_s, real_t omega1_s, 
                             real_t omega2_s) -> real_t;
  
  static auto dx_rsqr(real_t lambda, real_t press, 
                      real_t mbyr3) -> real_t;
  
  static auto dx_lambda(real_t press, real_t rho_e, 
                        real_t mbyr3) -> real_t;
  auto drsqr_omega1(real_t rsqr, real_t omega2) const -> real_t;
  
  static auto drsqr_omega2(real_t rsqr, real_t lambda, real_t rho, 
                           real_t hm1, real_t drsqr_omega1, 
                           real_t omega1) -> real_t;
  
  
  static auto ebnd_by_r3(real_t ybnd, real_t rsqr, real_t rho, 
                         real_t eps) -> real_t;

  static auto drsqr_ybnd(real_t ybnd, real_t rsqr, real_t lambda, 
                         real_t rho, real_t eps) -> real_t;

  static auto vol_by_r3(real_t yvol, real_t rsqr) -> real_t;
  
  static auto drsqr_yvol(real_t yvol, real_t rsqr, 
                         real_t lambda) -> real_t;


  void operator()(const state_t &s , state_t &dsdx, 
                       const real_t x) const;

  
  auto initial_data() const -> state_t;
  
  auto star(const state_t& surf) const -> spherical_star_info; 
};


}


#endif

