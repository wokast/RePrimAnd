#include "datastore.h"
#include "eos_barotr_file_impl.h"
#include "eos_barotr_table.h"
#include "eos_barotr_table_impl.h"

namespace EOS_Toolkit {
namespace implementations {

const std::string eos_barotr_table::datastore_id{ "barotr_table" };

struct reader_eos_barotr_table : reader_eos_barotr 
{
  eos_barotr load(const datasource g, const units& u) const final;
};

const bool eos_barotr_table::file_handler_registered { 
  registry_reader_eos_barotr::add(eos_barotr_table::datastore_id, 
                                  new reader_eos_barotr_table())
};

eos_barotr reader_eos_barotr_table::load(const datasource g, 
                                         const units& u) const
{
  
  bool isentropic = g["isentropic"];
  real_t poly_n   = g["poly_n"];

  std::vector<real_t> v_temp;
  if (g.has_data("temp")) {
    v_temp = g["temp"];
  }
  
  std::vector<real_t> v_ye;
  if (g.has_data("efr")) { 
    v_ye = g["efr"];
  }
  
  std::vector<real_t> v_rho = g["rmd"];
  std::vector<real_t> v_gm1 = g["gm1"];
  std::vector<real_t> v_eps = g["sed"];
  std::vector<real_t> v_p   = g["press"];
  std::vector<real_t> v_cs  = g["csnd"];

  std::size_t sz = v_rho.size();
  
  if ((v_gm1.size() != sz) || (v_eps.size() != sz) ||
      (v_p.size() != sz) || (v_cs.size() != sz) ||
      (!v_temp.empty() && v_temp.size() != sz) ||
      (!v_ye.empty() && v_ye.size() != sz))
  {
    throw std::runtime_error("Corrupt tabulated barotropic EOS file "
                             "(mismatching table sizes)"); 
  }

  std::vector<real_t> v_pbr(sz);
  std::vector<real_t> v_cs2(sz);

  for (std::size_t i=0; i < sz; ++i) {
    v_rho[i] /= u.density();
    v_p[i]   /= u.pressure();
    v_cs[i]  /= u.velocity();

    v_pbr[i] = v_p[i] / v_rho[i]; 
    v_cs2[i] = pow(v_cs[i], 2); 
  }

  
  return make_eos_barotr_table(v_gm1, v_rho, v_eps, v_pbr, v_cs2, 
                               v_temp, v_ye, isentropic, poly_n, u);
  
}
  
  
  
} 
}
