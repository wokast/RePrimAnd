#include <stdexcept>
#include <boost/math/tools/roots.hpp>
#include "spherical_stars_internals.h"


namespace EOS_Toolkit {
namespace details {
  
auto find_bulk_props(const spherical_star_profile& prf, 
                         real_t acc, std::size_t max_it) -> spherical_star_bulk
{
  auto froot = [&] (real_t rc) {
    real_t rho{ prf.state_from_rc(rc).rho() }; 
    real_t rho_avg{
      (rc > 0) ? (prf.mbary_from_rc(rc) / prf.pvol_from_rc(rc)) : rho
    };
    return rho_avg - 3. * rho;
  };
    
  real_t rmax{ prf.surf_circ_radius() };
  
  auto stopif = [&] (real_t ra, real_t rb) {
    return std::fabs(ra-rb) < rmax * acc;
  };
  
  boost::uintmax_t iters{ max_it };
  
  const auto rblk {
    boost::math::tools::toms748_solve(froot, 0., rmax, stopif, iters)
  };

  if (iters == max_it) {
    throw std::runtime_error("Root finding for bulk radius failed.");
  }

  const real_t blk_r{ (rblk.first + rblk.second) / 2. };
  const real_t blk_rho{ prf.state_from_rc(blk_r).rho() };
  const real_t blk_pv{ prf.pvol_from_rc(blk_r) };
  const real_t blk_mb{ prf.mbary_from_rc(blk_r) };
  
  return {blk_r, blk_rho, blk_pv, blk_mb}; 
}


}
}

