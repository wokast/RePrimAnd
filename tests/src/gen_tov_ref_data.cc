#include <boost/json.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>
#include <fstream>
#include <iomanip>
#include <array>
#include <stdexcept>
#include "config.h"
#include "test_config.h"
#include "unitconv.h"
#include "eos_barotropic.h"
#include "eos_barotr_file.h"
#include "eos_barotr_poly.h"
#include "spherical_stars.h"

using namespace std;
using namespace EOS_Toolkit;
using namespace boost::json;


void generate_ref_tov(eos_barotr eos, const real_t rho_cen, 
                      std::string tov_path, const units u)
{
  const real_t acc{ 1e-10 };
  const tov_acc_precise accs{acc, acc, acc, acc, 5000};
  
  const auto tov{ make_tov_star(eos, rho_cen, accs, true) };
  
  value props = {
    { "accuracy", acc },
    { "rho_cen", rho_cen },
    { "rho_cen_SI", rho_cen * u.density() },
    { "circ_radius", tov.circ_radius() },
    { "circ_radius_SI", tov.circ_radius() * u.length()},
    { "grav_mass", tov.grav_mass() },
    { "grav_mass_SI", tov.grav_mass() * u.mass()},
    { "bary_mass", tov.bary_mass() },
    { "bary_mass_SI", tov.bary_mass() * u.mass()},
    { "proper_volume", tov.proper_volume() },
    { "proper_volume_SI", tov.proper_volume() * u.volume()},
    { "moment_inertia", tov.moment_inertia() },
    { "moment_inertia_SI", tov.moment_inertia() * u.mom_inertia()},
    { "tidal_lambda", tov.deformability().lambda },
    { "tidal_k2", tov.deformability().k2 },
    { "bulk_radius", tov.bulk().circ_radius },
    { "bulk_radius_SI", tov.bulk().circ_radius * u.length()},
    { "bulk_bary_mass", tov.bulk().bary_mass },
    { "bulk_bary_mass_SI", tov.bulk().bary_mass * u.mass()},
    { "bulk_proper_volume", tov.bulk().proper_volume },
    { "bulk_proper_volume_SI", tov.bulk().proper_volume * u.volume() }
  };
  
  std::ofstream of(tov_path, ofstream::out);
  
  of << serialize(props) << endl;
  
}

int main() {
  auto u = units::geom_solar();
  const real_t ref_mass = 1.4;
  std::string eoslist[] = {"H4_Read_PP", "WFF1_Read_PP", 
                           "APR4_Read_PP", "MPA1_Read_PP",
                           "MS1_Read_PP"};
  for (auto s : eoslist) 
  {
    std::string eos_path{ std::string(PATH_TOV_EOS) + "/" 
                           + s + ".eos.h5" };
    std::string tov_path14{ std::string(PATH_TOV_REF) + "/ref_tov_m14_" 
                           + s + ".json" };
    std::string tov_pathmm{ std::string(PATH_TOV_REF) + "/ref_tov_max_" 
                           + s + ".json" };
    
    
    eos_barotr eos = load_eos_barotr(eos_path, u);
    
    const real_t rhomm0{ 1e17 / u.density()};
    const real_t rhomm1{ 1e19 / u.density()};
    const real_t rhomm{ 
      find_rhoc_tov_max_mass(eos, rhomm0, rhomm1) 
    };
    
    const real_t rho_ref{ 
      find_rhoc_tov_of_mass(eos, ref_mass, rhomm0, rhomm)
    };
    
    generate_ref_tov(eos, rho_ref, tov_path14, u);
    generate_ref_tov(eos, rhomm, tov_pathmm, u);
    
  }
}

