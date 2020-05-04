#include "hdf5imple.h"
#include "eos_barotr_file_impl.h"
#include "eos_barotr_poly.h"

namespace EOS_Toolkit {
namespace implementations {


struct reader_eos_barotr_poly : reader_eos_barotr 
{
  eos_barotr load(const h5grp& g, const units& u) const final;
};

const bool register_reader_eos_barotr_poly { 
  registry_reader_eos_barotr::add("barotr_poly", 
                                  new reader_eos_barotr_poly())
};

eos_barotr reader_eos_barotr_poly::load(const h5grp& g, 
                                        const units& u) const
{
  
  real_t poly_n  = get_attribute<real_t>(g, "poly_n");
  real_t rho_p   = get_attribute<real_t>(g, "rho_poly") / u.density();
  real_t rho_max = get_attribute<real_t>(g, "rho_max") / u.density();

  return make_eos_barotr_poly(poly_n, rho_p, rho_max);
}
  
  
  
} 
}
