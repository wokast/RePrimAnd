#define BOOST_TEST_MODULE TOV
#include <boost/test/unit_test.hpp>

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>

#include "test_utils.h"
#include <boost/format.hpp>
#include <boost/json.hpp>
#include "test_config.h"

#include "unitconv.h"
#include "eos_barotropic.h"
#include "eos_barotr_poly.h"
#include "eos_barotr_file.h"
#include "spherical_stars.h"


using boost::format;

using namespace std;
using namespace EOS_Toolkit;


auto load_json(std::string path) -> boost::json::value
{
  std::ifstream f(path);
  std::string s;
  if (f) {
    std::ostringstream ss;
    ss << f.rdbuf();
    s = ss.str();
  }
  return boost::json::parse(s);
}

auto check_tov_close(const boost::json::object& ref, 
                     const spherical_star& tov, 
                     real_t tol_tov, real_t tol_def) -> bool
{
  failcount hope("TOV solutions agree within tolerance");
  
  hope.isclose(tov.grav_mass(), ref.at("grav_mass").get_double(), 
               tol_tov, 0, "grav. mass");

  hope.isclose(tov.bary_mass(), ref.at("bary_mass").get_double(), 
               tol_tov, 0, "bary. mass");
  hope.isclose(tov.circ_radius(), ref.at("circ_radius").get_double(), 
               tol_tov, 0, "circ. radius");
  hope.isclose(tov.proper_volume(), ref.at("proper_volume").get_double(), 
               3. * tol_tov, 0, "proper volume");
  hope.isclose(tov.moment_inertia(), 
               ref.at("moment_inertia").get_double(), 
               tol_tov, 0, "moment of inertia");
               
  hope.isclose(tov.deformability().lambda,  
               ref.at("tidal_lambda").get_double(), 
               tol_def, 0, "tidal deformability");
  hope.isclose(tov.deformability().k2,  
               ref.at("tidal_k2").get_double(), 
               tol_def, 0, "love number");
               
  hope.isclose(tov.bulk().circ_radius, ref.at("bulk_radius").get_double(), 
               tol_tov, 0, "bulk radius");

  hope.isclose(tov.bulk().bary_mass, ref.at("bulk_bary_mass").get_double(), 
               4. * tol_tov, 0, "bulk baryonic mass");
               
  hope.isclose(tov.bulk().proper_volume, 
               ref.at("bulk_proper_volume").get_double(), 
               3. * tol_tov, 0, "bulk proper volume");

  return hope;
}



bool check_ref_tov(std::string eospath, std::string jsonpath)
{

  auto ref{ load_json(jsonpath).get_object() };

  auto uc{ units::geom_solar() };
  
  eos_barotr eos{ load_eos_barotr(eospath, uc) };
  const real_t rho_cen{ ref["rho_cen"].get_double() };


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
                           + s + ".json" };
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
