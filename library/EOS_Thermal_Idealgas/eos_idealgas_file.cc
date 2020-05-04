#include "hdf5imple.h"
#include "eos_thermal_file_impl.h"
#include "eos_idealgas.h"

namespace EOS_Toolkit {
namespace implementations {


struct reader_eos_thermal_idealgas : reader_eos_thermal 
{
  eos_thermal load(const h5grp& g, const units& u) const final;
};

const bool register_reader_eos_thermal_idealgas { 
  registry_reader_eos_thermal::add("thermal_idealgas", 
                                  new reader_eos_thermal_idealgas())
};

eos_thermal reader_eos_thermal_idealgas::load(const h5grp& g, 
                                           const units& u) const
{
  real_t n_adiab   = get_attribute<real_t>(g, "adiab_index");
  real_t eps_max   = get_attribute<real_t>(g, "eps_max");
  real_t rho_max   = get_attribute<real_t>(g, "rho_max") / u.density();
                    
  return make_eos_idealgas(n_adiab, eps_max, rho_max);
}
  
  
  
} //namespace implementations
} //namespace EOS_Toolkit
