#include "datastore.h"
#include "eos_barotr_file_impl.h"
#include "eos_barotr_pwpoly.h"
#include "eos_barotr_pwpoly_impl.h"

namespace EOS_Toolkit {
namespace implementations {

const std::string eos_barotr_pwpoly::datastore_id{ "barotr_pwpoly" };

struct reader_eos_barotr_pwpoly : reader_eos_barotr 
{
  eos_barotr load(const datasource g, const units& u) const final;
};

const bool eos_barotr_pwpoly::file_handler_registered {
  registry_reader_eos_barotr::add(eos_barotr_pwpoly::datastore_id, 
                                  new reader_eos_barotr_pwpoly())
};

eos_barotr reader_eos_barotr_pwpoly::load(const datasource g, 
                                          const units& u) const
{
  if (g.has_data("eos_type")) { //missing in old format
    std::string id = g["eos_type"];
    if (id != eos_barotr_pwpoly::datastore_id)
    {
      throw std::runtime_error("eos_barotr_pwpoly: trying to load from "
                          "stored EOS of different type");
    }
  }  
  
  real_t rho_p_si   = g["rho_poly"];
  real_t rho_p      = rho_p_si / u.density();
  
  real_t rho_max_si = g["rho_max"];
  real_t rho_max    = rho_max_si / u.density();
  
  std::vector<real_t> v_rho_b = g["rho_bound"];
  std::vector<real_t> v_gamma = g["gamma"];
  
  for (auto &r : v_rho_b) r /= u.density();
  
  return make_eos_barotr_pwpoly(rho_p, v_rho_b, v_gamma, rho_max, u); 
} 


void eos_barotr_pwpoly::save(datasink g) const
{
  auto u = units_to_SI();
  g["eos_type"] = datastore_id;
  
  g["rho_poly"] = segments[0].rmd_p * u.density();  
  g["rho_max"]  = range_rho().max() * u.density();
  
  std::vector<real_t> v_rho_b_si;
  std::vector<real_t> v_gamma;
  for (auto s : segments) {
    v_rho_b_si.push_back( s.rmd0 * u.density() );
    v_gamma.push_back( s.gamma );
  }
  g["rho_bound"] = v_rho_b_si;
  g["gamma"]     = v_gamma;
}  
  

} //namespace implementations 
} // namespace EOS_Toolkit 
