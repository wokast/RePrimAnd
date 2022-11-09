#ifndef TOV_SEQS_IMPL_H
#define TOV_SEQS_IMPL_H
#include <memory>
#include <vector>
#include <cassert>

#include "star_sequence.h"
#include "interpol.h"
#include "unitconv.h"
#include "datastore.h"

namespace EOS_Toolkit {



namespace detail {


class star_seq_impl {
  public:

  using spline_t = interpolator;  
  using range_t  = interval<real_t>;
  
  star_seq_impl(spline_t mg_gm1_, spline_t mb_gm1_, 
                spline_t rc_gm1_, spline_t mi_gm1_, 
                spline_t lt_gm1_, range_t rg_gm1_,
                units u_);

  static auto from_vector(std::vector<real_t> mg, 
     std::vector<real_t> mb, std::vector<real_t> rc, 
     std::vector<real_t> mi, std::vector<real_t> lt, 
     range_t rg_gm1, units u)
  -> std::shared_ptr<star_seq_impl>;
  
  auto grav_mass_from_center_gm1(real_t gm1c) const -> real_t; 
  auto bary_mass_from_center_gm1(real_t gm1c) const -> real_t; 
  auto circ_radius_from_center_gm1(real_t gm1c) const -> real_t; 
  auto moment_inertia_from_center_gm1(real_t gm1c) const -> real_t; 
  auto lambda_tidal_from_center_gm1(real_t gm1c) const -> real_t; 
  
  auto range_center_gm1() const -> range_t;
  auto contains_gm1(real_t gm1c) const -> bool;
  
  void save(datasink s) const;
  
  auto static from_datasource(datasource g, units u) 
  -> std::shared_ptr<star_seq_impl>;
  
  auto units_to_SI() const -> const units&;

  private:
  
  spline_t mg_gm1;
  spline_t mb_gm1;
  spline_t rc_gm1;
  spline_t mi_gm1;
  spline_t lt_gm1;
  
  const range_t rg_gm1;
  
  units u;
};

class star_branch_impl {
  public:
  
  using spline_t = interpolator;  
  using range_t  = interval<real_t>;
  
  star_branch_impl(range_t rg_gm1_, 
                   spline_t xg_mg_, real_t gm1_ref_,
                   bool incl_max_, units u_);
  
  
  auto center_gm1_from_grav_mass(real_t mg) const -> real_t; 
  
  auto range_center_gm1() const -> range_t;
  auto contains_gm1(real_t gm1c) const -> bool;

  auto range_grav_mass() const -> range_t;
  auto contains_grav_mass(real_t mg) const -> bool;
  
  auto includes_maximum() const -> bool;
  auto grav_mass_maximum() const -> real_t;
  //~ auto bary_mass_maximum() const -> real_t;
  auto center_gm1_maximum() const -> real_t;

  void save(datasink s) const;

  auto static from_datasource(datasource g, units u) 
  -> std::shared_ptr<star_branch_impl>;

  auto static xg_from_gm1(real_t gm1, real_t gm1_ref) -> real_t;
  private:
  
  auto gm1_from_xg(real_t xg) const -> real_t;

  range_t rg_gm1;
  spline_t xg_mg;  
  const real_t gm1_ref;
  const bool incl_max;
  units u;
};


}
}

#endif
