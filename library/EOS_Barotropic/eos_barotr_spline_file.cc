#include <boost/optional.hpp>
#include "datastore.h"
#include "interpol_logspl.h"
#include "eos_barotr_file_impl.h"
#include "eos_barotr_spline_impl.h"

namespace EOS_Toolkit {
namespace implementations {

const std::string eos_barotr_spline::datastore_id{ "barotr_spline_v2" };

struct reader_eos_barotr_spline : reader_eos_barotr 
{
  eos_barotr load(const datasource g, const units& u) const final;
};

const bool eos_barotr_spline::file_handler_registered { 
  registry_reader_eos_barotr::add(eos_barotr_spline::datastore_id, 
                                  new reader_eos_barotr_spline())
};

eos_barotr reader_eos_barotr_spline::load(const datasource g, 
                                          const units& u) const
{
  using lgspl_t = detail::interpol_logspl_impl;
  using lglgspl_t = detail::interpol_llogspl_impl;
  using opt_t = boost::optional<lgspl_t>;
  
  std::string id = g["eos_type"];
  if (id != eos_barotr_spline::datastore_id)
  {
    throw std::runtime_error("eos_barotr_spline: trying to load from "
                             "stored EOS of different type");
  }
  
  bool isentropic = g["isentropic"];
  
  auto poly = eos_barotr_gpoly::from_datasource(g / "eos_gpoly", u);

  lglgspl_t sgm1_si   = g["gm1_from_rho"];
  lglgspl_t srho_si   = g["rho_from_gm1"];
  lgspl_t seps        = g["eps_from_gm1"];
  lgspl_t shm1        = g["hm1_from_gm1"];
  lglgspl_t spress_si = g["press_from_gm1"];
  lgspl_t scsnd_si    = g["csnd_from_rho"];
  
  opt_t stemp_mev     = g["temp_from_gm1"];
  opt_t sefrac        = g["efrac_from_gm1"];
  
  lglgspl_t sgm1   = sgm1_si.rescale_x(1./u.density());
  lglgspl_t srho   = srho_si / u.density();
  lglgspl_t spress = spress_si / u.pressure();
  lgspl_t scsnd    = scsnd_si.rescale_x(1./u.density()) / u.velocity();

  return eos_barotr{ 
    std::make_shared<eos_barotr_spline>(sgm1, srho, seps, spress, 
                   shm1, scsnd, stemp_mev, sefrac, isentropic, poly) 
  };
}


void eos_barotr_spline::save(datasink g) const
{
  auto u = units_to_SI();
  g["eos_type"] = datastore_id;
  
  poly.save(g / "eos_gpoly");
  g["isentropic"] = isentropic;
  g["gm1_from_rho"] = gm1_rho.rescale_x(u.density());
  g["rho_from_gm1"] = rho_gm1 * u.density();
  g["eps_from_gm1"] = eps_gm1;
  g["hm1_from_gm1"] = hm1_gm1;
  g["press_from_gm1"] = p_gm1 * u.pressure();
  g["csnd_from_rho"] = csnd_rho.rescale_x(u.density()) * u.velocity();
  
  
  if (!zerotemp) 
  {
    g["temp_from_gm1"] = *temp_gm1;
  }
  
  if (efrac_gm1)
  {
    g["efrac_from_gm1"] = *efrac_gm1;
  }
  
}  
  
  
} 
}
