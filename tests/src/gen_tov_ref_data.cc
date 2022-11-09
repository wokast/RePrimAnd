#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>
#include <fstream>
#include <iomanip>
#include <array>
#include <stdexcept>
#include "config.h"
#include "test_config.h"
#include "unitconv.h"
#include "hdf5store.h"
#include "eos_barotropic.h"
#include "eos_barotr_file.h"
#include "eos_barotr_poly.h"
#include "spherical_stars.h"
#include "star_sequence.h"
#include "star_seq_file.h"

using namespace std;
using namespace EOS_Toolkit;


void generate_ref_tov(eos_barotr eos, const real_t rho_cen, 
                      std::string tov_path, const units u)
{
  const real_t acc{ 1e-10 };
  const tov_acc_precise accs{acc, acc, acc, acc, 5000};
  
  const auto tov{ make_tov_star(eos, rho_cen, accs, true) };
  
  auto s = make_hdf5_file_sink(tov_path);
    
  s["accuracy"] = acc;
  s["rho_cen"] = rho_cen;
  s["rho_cen_SI"] = rho_cen * u.density();
  s["circ_radius"] = tov.circ_radius();
  s["circ_radius_SI"] = tov.circ_radius() * u.length();
  s["grav_mass"] = tov.grav_mass();
  s["grav_mass_SI"] = tov.grav_mass() * u.mass();
  s["bary_mass"] = tov.bary_mass();
  s["bary_mass_SI"] = tov.bary_mass() * u.mass();
  s["proper_volume"] = tov.proper_volume();
  s["proper_volume_SI"] = tov.proper_volume() * u.volume();
  s["moment_inertia"] = tov.moment_inertia();
  s["moment_inertia_SI"] = tov.moment_inertia() * u.mom_inertia();
  s["tidal_lambda"] = tov.deformability().lambda;
  s["tidal_k2"] = tov.deformability().k2;
  s["bulk_radius"] = tov.bulk().circ_radius;
  s["bulk_radius_SI"] = tov.bulk().circ_radius * u.length();
  s["bulk_bary_mass"] = tov.bulk().bary_mass;
  s["bulk_bary_mass_SI"] = tov.bulk().bary_mass * u.mass();
  s["bulk_proper_volume"] = tov.bulk().proper_volume;
  s["bulk_proper_volume_SI"] = tov.bulk().proper_volume * u.volume();

}

auto get_eos_by_name(std::string s, units u)
->eos_barotr
{
    std::string eos_path{ std::string(PATH_TOV_EOS) + "/" 
                           + s + ".eos.h5" };
    return load_eos_barotr(eos_path, u);
}

void make_ref_tovs()
{
  auto u = units::geom_solar();
  const real_t ref_mass = 1.4;
  std::string eoslist[] = {"H4_Read_PP", "H4_Read_PP.spline", 
                           "WFF1_Read_PP", 
                           "APR4_Read_PP", "MPA1_Read_PP",
                           "MS1_Read_PP"};
  for (auto s : eoslist) 
  {
    std::string tov_path14{ std::string(PATH_TOV_REF) + "/ref_tov_m14_" 
                           + s + ".h5" };
    std::string tov_pathmm{ std::string(PATH_TOV_REF) + "/ref_tov_max_" 
                           + s + ".h5" };
  
    auto eos{ get_eos_by_name(s, u) };
    
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

void make_ref_seqs()
{
  auto u = units::geom_solar();
  std::string eoslist[] = {"H4_Read_PP", "H4_Read_PP.spline",
                           "WFF1_Read_PP", 
                           "APR4_Read_PP", "MPA1_Read_PP",
                           "MS1_Read_PP", "APR4_EPP"};
  for (auto s : eoslist) 
  {
    std::string seq_path{ std::string(PATH_TOV_REF) + "/ref_tovseq_" 
                           + s + ".h5" };
    
    auto eos{ get_eos_by_name(s, u) };

    const tov_acc_simple acc{1e-8, 1e-6};

    auto seq = make_tov_branch_stable(eos, acc);
    
    save_star_branch(seq_path, seq);
  }
}

int main() {
  make_ref_tovs();
  make_ref_seqs();
}

