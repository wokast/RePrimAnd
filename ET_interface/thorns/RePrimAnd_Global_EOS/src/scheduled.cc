#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "reprimand_global_eos.h"
#include "reprimand/unitconv.h"
#include "reprimand/eos_idealgas.h"
#include "reprimand/eos_hybrid.h"
#include "reprimand/eos_thermal_file.h"
#include "reprimand/eos_barotr_poly.h"
#include "reprimand/eos_barotr_file.h"


using namespace EOS_Toolkit;
using namespace RePrimAnd_Global_EOS;

extern "C" int RePrimAnd_Evol_EOS_Register_IdealGas()
{
  DECLARE_CCTK_PARAMETERS;

  try { 
    auto eos = make_eos_idealgas(idealgas_n, eps_max, rho_max);
    global_eos_thermal::set_eos(eos);
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}

extern "C" int RePrimAnd_Evol_EOS_Register_Hybrid()
{
  DECLARE_CCTK_PARAMETERS;

  try { 
    auto eos_cold = global_eos_cold::get_eos();
    auto eos = make_eos_hybrid(eos_cold, hybrid_gamma_th, 
                               eps_max, rho_max);
    global_eos_thermal::set_eos(eos);
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}

extern "C" int RePrimAnd_Evol_EOS_Register_FromFile()
{
  DECLARE_CCTK_PARAMETERS;

  try { 
    auto eos_units = units::geom_ulength(geom_units_length);
    auto eos = load_eos_thermal(thermal_eos_file, eos_units);
    global_eos_thermal::set_eos(eos);
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}

extern "C" int RePrimAnd_Initial_EOS_Register_Polytrope()
{
  DECLARE_CCTK_PARAMETERS;

  try { 
    auto u = units::geom_ulength(geom_units_length);
    auto eos = make_eos_barotr_poly(poly_n, poly_rho / u.density(), 
                                    rho_max);
    global_eos_cold::set_eos(eos);
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}

extern "C" int RePrimAnd_Initial_EOS_Register_FromFile()
{
  DECLARE_CCTK_PARAMETERS;

  try { 
    auto eos_units = units::geom_ulength(geom_units_length);
    auto eos = load_eos_barotr(cold_eos_file, eos_units);
    global_eos_cold::set_eos(eos);
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}


