#include "datastore.h"
#include "eos_thermal_file_impl.h"
#include "eos_idealgas.h"
#include "eos_idealgas_impl.h"

namespace EOS_Toolkit {
namespace implementations {

const std::string eos_idealgas::datastore_id{"thermal_idealgas"};

struct reader_eos_thermal_idealgas : reader_eos_thermal 
{
  eos_thermal load(const datasource g, const units& u) const final;
};

const bool eos_idealgas::file_handler_registered { 
  registry_reader_eos_thermal::add(eos_idealgas::datastore_id, 
                                  new reader_eos_thermal_idealgas())
};

eos_thermal reader_eos_thermal_idealgas::load(const datasource g, 
                                           const units& u) const
{
  real_t n_adiab    = g["adiab_index"];
  real_t eps_max    = g["eps_max"];
  
  real_t rho_max_si = g["rho_max"];
  real_t rho_max    = rho_max_si / u.density();
                    
  return make_eos_idealgas(n_adiab, eps_max, rho_max, u);
}
  

void eos_idealgas::save(datasink g) const
{
  auto u = units_to_SI();
  g["eos_type"] = datastore_id;
  
  g["adiab_index"] = n_index;
  g["eps_max"] = rgeps.max();
  g["rho_max"] = rgrho.max() *  u.density();
}    
  
} //namespace implementations
} //namespace EOS_Toolkit
