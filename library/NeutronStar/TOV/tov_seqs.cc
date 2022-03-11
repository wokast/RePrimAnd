#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>
#include <stdexcept>
#include "spherical_stars.h"


using namespace std;


namespace EOS_Toolkit {

auto find_rhoc_tov_max_mass(eos_barotr eos, 
                       const real_t rhobr0, const real_t rhobr1,
                       const int bits, const real_t acc, 
                       unsigned int max_steps)
-> real_t
{
  const real_t rho0 { eos.range_rho().limit_to(rhobr0) };
  const real_t rho1 { eos.range_rho().limit_to(rhobr1) };
  
  tov_acc_simple accs{acc};
  
  auto fmin = [&] (real_t rho_c) {
    const auto tov{ 
      get_tov_star_properties(eos, rho_c, accs, false, false) 
    };
    return -tov.grav_mass();
  };
  
  boost::uintmax_t iters{ max_steps };
  const auto r{ 
    boost::math::tools::brent_find_minima(fmin, rho0, rho1, 
                                          bits, iters)
  };
  
  if (iters >= max_steps)
  {
    throw std::runtime_error("TOV maximum mass not found");
  }
  
  return r.first;
  
}



auto find_rhoc_tov_of_mass(eos_barotr eos, real_t mg, 
                        const real_t rhobr0, const real_t rhobr1,
                        real_t acc, unsigned int max_steps)
-> real_t


{
  tov_acc_simple accs{acc};
  
  auto froot = [&] (real_t rho_c) {
    auto tov = get_tov_star_properties(eos, rho_c, accs, false, false);
    return tov.grav_mass() - mg;
  };

  auto stopif = [&] (real_t l, real_t r) {
    return fabs(l-r) < fabs(l+r) * acc;
  };
  
  boost::uintmax_t iters{ max_steps };
    
  auto res = boost::math::tools::toms748_solve(froot, rhobr0, rhobr1, 
                                               stopif, iters);
  
  if (iters >= max_steps)
  {
    throw std::runtime_error("TOV model with mass: root finding failed");
  }
  
  return res.first;
}


}


