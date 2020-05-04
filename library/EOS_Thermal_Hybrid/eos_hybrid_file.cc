#include "hdf5imple.h"
#include "eos_barotr_file.h"
#include "eos_thermal_file_impl.h"
#include "eos_hybrid.h"

namespace EOS_Toolkit {
namespace implementations {


struct reader_eos_thermal_hybrid : reader_eos_thermal 
{
  eos_thermal load(const h5grp& g, const units& u) const final;
};

const bool register_reader_eos_thermal_hybrid { 
  registry_reader_eos_thermal::add("thermal_hybrid", 
                                  new reader_eos_thermal_hybrid())
};

eos_thermal reader_eos_thermal_hybrid::load(const h5grp& g, 
                                           const units& u) const
{
  real_t gamma_th  = get_attribute<real_t>(g, "gamma_th");
  real_t eps_max   = get_attribute<real_t>(g, "eps_max");

  h5grp g2(g, "eos_cold");
  
  auto eos_cold = ::EOS_Toolkit::detail::load_eos_barotr(g2, u);

  real_t rho_max = eos_cold.range_rho().max();
                                
  return make_eos_hybrid(eos_cold, gamma_th, eps_max, rho_max);
}
  
  
  
} 
}
