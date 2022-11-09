#ifndef SPHERICAL_STARS_INTERNALS_H
#define SPHERICAL_STARS_INTERNALS_H

#include "spherical_stars.h"
#include "tov_ode.h"
#include "interpol.h"


namespace EOS_Toolkit {
namespace details {
  
  
class tov_profile : public spherical_star_profile {
  const interpolator lambda_rsqr;
  const interpolator delta_nu_rsqr;
  const interpolator ybnd_rsqr;
  const interpolator yvol_rsqr;
  const real_t gm1_c;
  const real_t nu_c;
  const real_t mgrav;
  const real_t ebind;

  void validate_rc(real_t rc) const;

  auto nu_from_rc_outside(real_t rc) const -> real_t;
  auto pvol_vacuum(real_t rc) const -> real_t;

  public:
  using vec_t = std::vector<real_t>;

  tov_profile() = delete;
  tov_profile(tov_profile&&) = default;
  tov_profile(const tov_profile&) = default;
  
  tov_profile& operator=(tov_profile&&) = default;
  tov_profile& operator=(const tov_profile&) = default;
  
  tov_profile(eos_barotr eos_, const spherical_star_info &p_,
              vec_t rsqr_, vec_t delta_nu_, vec_t lambda_,
              vec_t ybnd_, vec_t yvol_);
  
  auto center_gm1() const -> real_t override;
  auto nu_from_rc(real_t rc) const -> real_t override;
  auto lambda_from_rc(real_t rc) const -> real_t override;
  auto gm1_from_rc(real_t rc) const -> real_t override;
  auto mbary_from_rc(real_t rc) const -> real_t override;
  auto pvol_from_rc(real_t rc) const -> real_t override;
};


auto find_bulk_props(const spherical_star_profile& prf, real_t acc, 
               std::size_t max_it=30) -> spherical_star_bulk;
                         
}
}




#endif
