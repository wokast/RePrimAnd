#include "hdf5imple.h"
#include "eos_barotr_file_impl.h"
#include "eos_barotr_table.h"

namespace EOS_Toolkit {
namespace implementations {


struct reader_eos_barotr_table : reader_eos_barotr 
{
  eos_barotr load(const h5grp& g, const units& u) const final;
};

const bool register_reader_eos_barotr_table { 
  registry_reader_eos_barotr::add("barotr_table", 
                                  new reader_eos_barotr_table())
};

eos_barotr reader_eos_barotr_table::load(const h5grp& g, 
                                         const units& u) const
{
  
  bool isentropic = get_attribute<bool>(g, "isentropic");
  real_t poly_n   = get_attribute<real_t>(g, "poly_n");

  std::vector<real_t> v_temp;
  h5dset d_temp(g, "temp");
  if (d_temp) {
    read(d_temp, v_temp);
  }
  std::vector<real_t> v_ye;
  h5dset d_ye(g, "efr");
  if (d_ye) { 
    read(d_ye, v_ye);
  }
  
  auto v_rho = get_vector<real_t>(g, "rmd");
  auto v_gm1 = get_vector<real_t>(g, "gm1");
  auto v_eps = get_vector<real_t>(g, "sed");
  auto v_p   = get_vector<real_t>(g, "press");
  auto v_cs  = get_vector<real_t>(g, "csnd");

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
                               v_temp, v_ye, isentropic, poly_n);
  
}
  
  
  
} 
}
