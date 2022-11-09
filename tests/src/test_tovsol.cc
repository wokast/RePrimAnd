#define BOOST_TEST_MODULE TOV
#include <boost/test/unit_test.hpp>

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>

#include "test_utils.h"
#include "test_config.h"

#include "unitconv.h"
#include "hdf5store.h"
#include "eos_barotropic.h"
#include "eos_barotr_poly.h"
#include "eos_barotr_file.h"
#include "spherical_stars.h"



using namespace std;
using namespace EOS_Toolkit;


auto check_tov_close(datasource s, 
                     const spherical_star& tov, 
                     real_t tol_tov, real_t tol_def) -> bool
{
  failcount hope("TOV solutions agree within tolerance");
  
  hope.isclose(tov.grav_mass(), real_t(s["grav_mass"]), 
               tol_tov, 0, "grav. mass");

  hope.isclose(tov.bary_mass(), real_t(s["bary_mass"]), 
               tol_tov, 0, "bary. mass");
  hope.isclose(tov.circ_radius(), real_t(s["circ_radius"]), 
               tol_tov, 0, "circ. radius");
  hope.isclose(tov.proper_volume(), real_t(s["proper_volume"]), 
               3. * tol_tov, 0, "proper volume");
  hope.isclose(tov.moment_inertia(), 
               real_t(s["moment_inertia"]), 
               tol_tov, 0, "moment of inertia");
               
  hope.isclose(tov.deformability().lambda,  
               real_t(s["tidal_lambda"]), 
               tol_def, 0, "tidal deformability");
  hope.isclose(tov.deformability().k2,  
               real_t(s["tidal_k2"]), 
               tol_def, 0, "love number");
               
  hope.isclose(tov.bulk().circ_radius, real_t(s["bulk_radius"]), 
               tol_tov, 0, "bulk radius");

  hope.isclose(tov.bulk().bary_mass, real_t(s["bulk_bary_mass"]), 
               4. * tol_tov, 0, "bulk baryonic mass");
               
  hope.isclose(tov.bulk().proper_volume, 
               real_t(s["bulk_proper_volume"]), 
               3. * tol_tov, 0, "bulk proper volume");

  return hope;
}



bool check_ref_tov(std::string eospath, std::string refpath)
{

  auto ref{ make_hdf5_file_source(refpath) };

  auto uc{ units::geom_solar() };
  
  eos_barotr eos{ load_eos_barotr(eospath, uc) };
  const real_t rho_cen{ real_t(ref["rho_cen"]) };


  const tov_acc_simple accs{1e-8, 1e-6, 500}; 
  const real_t tol_tov{ 1e2 * accs.tov};
  const real_t tol_def{ 5e-5 }; // limited by nsamp=500
  
  auto tov{ make_tov_star(eos, rho_cen, accs, true) }; 

  return check_tov_close(ref, tov, tol_tov, tol_def);

}


BOOST_AUTO_TEST_CASE( test_tovsol_tabeos )
{
  failcount hope("TOV solution for one tabulated EOS yields "
                 "known expected values for one density.");

  std::string eoslist[] = {
    "H4_Read_PP", 
    "H4_Read_PP.spline", 
    "WFF1_Read_PP", 
    "APR4_Read_PP", 
    "MPA1_Read_PP",
    "MS1_Read_PP"
    };
  for (auto s : eoslist) 
  {
    std::string eos_path{ std::string(PATH_TOV_EOS) + "/" 
                           + s + ".eos.h5" };
    std::string tov_path14{ std::string(PATH_TOV_REF) + "/ref_tov_m14_" 
                           + s + ".h5" };
    if (!hope(check_ref_tov(eos_path, tov_path14), 
       "Good accuracy solution within tolerance to reference")) 
    {
      hope.postmortem(tov_path14);
    }
  }
}
    

BOOST_AUTO_TEST_CASE( test_tovsol_basic )
{
  failcount hope("Finding a TOV solution works at all.");

  auto u = units::geom_solar();

  const real_t n_poly      = 1;
  const real_t rmd_poly    = 6.176e+18 / u.density();
  const real_t rhomax_poly = 1E40 / u.density();
  auto eos = make_eos_barotr_poly(n_poly, rmd_poly, rhomax_poly);

  const real_t tov_cnt_rmd       = 7.9056e+17 / u.density();
  
  const tov_acc_simple accs{1e-10, 1e-10, 500};
  
  hope.nothrow("Find TOV solution", 
               [&] () {make_tov_star(eos, tov_cnt_rmd, accs, true, true);}
              );
  
  hope.nothrow("Find TOV solution", 
               [&] () {make_tov_star(eos, tov_cnt_rmd, accs, false, false);}
              );

  auto tov = get_tov_star_properties(eos, tov_cnt_rmd, accs, true); 

  hope.istrue(tov.has_deform(), "tidal deformability available");
  hope.istrue(tov.has_bulk(), "bulk properties available");
}
