#include <string>
#include "datastore.h"
#include "eos_barotr_file_impl.h"
#include "eos_barotr_poly.h"
#include "eos_barotr_poly_impl.h"

namespace EOS_Toolkit {
namespace implementations {

const std::string eos_barotr_poly::datastore_id{ "barotr_poly" };

struct reader_eos_barotr_poly : reader_eos_barotr 
{
  eos_barotr load(const datasource g, const units& u) const final;
};

const bool eos_barotr_poly::file_handler_registered { 
  registry_reader_eos_barotr::add(eos_barotr_poly::datastore_id, 
                                  new reader_eos_barotr_poly())
};

eos_barotr reader_eos_barotr_poly::load(const datasource g, 
                                        const units& u) const
{
  if (g.has_data("eos_type")) { //missing in old format
    std::string id = g["eos_type"];
    if (id != eos_barotr_poly::datastore_id)
    {
      throw std::runtime_error("eos_barotr_poly: trying to load from "
                          "stored EOS of different type");
    }
  }  
  real_t poly_n  = g["poly_n"];
  
  real_t rho_p_si   = g["rho_poly"];
  real_t rho_p      = rho_p_si / u.density();
  
  real_t rho_max_si = g["rho_max"];
  real_t rho_max    = rho_max_si / u.density();

  return make_eos_barotr_poly(poly_n, rho_p, rho_max, u);
}


void eos_barotr_poly::save(datasink g) const
{
  auto u = units_to_SI();
  g["eos_type"] = datastore_id;
  
  g["poly_n"] = n;
  g["rho_poly"] = rmd_p * u.density();
  g["rho_max"] = range_rho().max() * u.density();
}  
  
  
} 
}
