#include "datastore.h"
#include "eos_barotr_file.h"
#include "eos_thermal_file_impl.h"
#include "eos_hybrid.h"
#include "eos_hybrid_impl.h"

namespace EOS_Toolkit {
namespace implementations {

const std::string eos_hybrid::datastore_id{"thermal_hybrid"};

struct reader_eos_thermal_hybrid : reader_eos_thermal 
{
  eos_thermal load(const datasource g, const units& u) const final;
};

const bool eos_hybrid::file_handler_registered { 
  registry_reader_eos_thermal::add(eos_hybrid::datastore_id, 
                                  new reader_eos_thermal_hybrid())
};

eos_thermal reader_eos_thermal_hybrid::load(const datasource g, 
                                           const units& u) const
{
  real_t gamma_th  = g["gamma_th"];
  real_t eps_max   = g["eps_max"];

  auto g2 = g / "eos_cold";
  
  auto eos_cold = ::EOS_Toolkit::detail::load_eos_barotr(g2, u);

  real_t rho_max = eos_cold.range_rho().max();
                                
  return make_eos_hybrid(eos_cold, gamma_th, eps_max, rho_max);
}
  
 
void eos_hybrid::save(datasink g) const
{
  g["eos_type"] = datastore_id;
  g["gamma_th"] = gamma_th;
  g["eps_max"] = eps_max;
  ::EOS_Toolkit::detail::save_eos_barotr(g / "eos_cold", eos_c);
}     
  
} 
}
