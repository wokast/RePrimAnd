#include "hdf5imple.h"
#include "eos_barotr_file_impl.h"
#include "eos_barotr_pwpoly.h"

namespace EOS_Toolkit {
namespace implementations {


struct reader_eos_barotr_pwpoly : reader_eos_barotr 
{
  eos_barotr load(const h5grp& g, const units& u) const final;
};

const bool register_reader_eos_barotr_pwpoly { 
  registry_reader_eos_barotr::add("barotr_pwpoly", 
                                  new reader_eos_barotr_pwpoly())
};

eos_barotr reader_eos_barotr_pwpoly::load(const h5grp& g, 
                                          const units& u) const
{
  
  real_t rho_p    = get_attribute<real_t>(g, "rho_poly") / u.density();
  real_t rho_max  = get_attribute<real_t>(g, "rho_max") / u.density();
  
  auto v_rho_b = get_vector<real_t>(g, "rho_bound");
  auto v_gamma = get_vector<real_t>(g, "gamma");
  
  for (auto &r : v_rho_b) r /= u.density();
  
  return make_eos_barotr_pwpoly(rho_p, v_rho_b, v_gamma, rho_max); 

} 


} //namespace implementations 
} // namespace EOS_Toolkit 
