#define BOOST_TEST_MODULE NSSEQS

#include <cstdio>
#include <boost/test/unit_test.hpp>
#include "test_utils.h"
#include <boost/format.hpp>

#include "test_config.h"

#include "unitconv.h"
#include "eos_barotropic.h"
#include "eos_barotr_file.h"
#include "star_sequence.h"
#include "star_seq_file.h"


using namespace std;
using namespace EOS_Toolkit;



bool compare_star_seqs(star_branch s1, star_branch s2, 
                       std::size_t nsamp,
                       real_t err_m_gm1, real_t err_rc_gm1,
                       real_t err_mi_gm1, real_t err_lt_gm1, 
                       real_t err_gm1_mg, real_t err_m_mg, 
                       real_t err_rc_mg, real_t err_mi_mg, 
                       real_t err_lt_mg, real_t err_rg)
{
  failcount hope("Check that two star sequences agree");
  
  hope(s1.includes_maximum() == s2.includes_maximum(),
       "seqs agree about max mass model is contained");
  
  hope.isclose(s1.range_center_gm1().max(), 
               s2.range_center_gm1().max(), err_rg, 0, 
               "Ranges g-1 upper");
  hope.isclose(s1.range_center_gm1().min(), 
               s2.range_center_gm1().min(), err_rg, 0, 
               "Ranges g-1 lower");
  hope.isclose(s1.range_grav_mass().max(), 
               s2.range_grav_mass().max(), err_rg, 0, 
               "Ranges mg upper");
  hope.isclose(s1.range_grav_mass().min(), 
               s2.range_grav_mass().min(), err_rg, 0, 
               "Ranges mg lower");
               
  auto rg_gm1{ intersect(s1.range_center_gm1(), 
                         s2.range_center_gm1()) };
  auto rg_mg{ intersect(s1.range_grav_mass(), 
                         s2.range_grav_mass()) };
  
  auto igm1 = linear_spacing(rg_gm1.min(), rg_gm1.max(), nsamp);
  for (real_t gm1 : igm1) 
  {
    hope.isclose(s1.grav_mass_from_center_gm1(gm1), 
                 s2.grav_mass_from_center_gm1(gm1), 
                 err_m_gm1, 0., "mg from g-1");
                 
    hope.isclose(s1.bary_mass_from_center_gm1(gm1), 
                 s2.bary_mass_from_center_gm1(gm1), 
                 err_m_gm1, 0., "mb from g-1");

    hope.isclose(s1.circ_radius_from_center_gm1(gm1), 
                 s2.circ_radius_from_center_gm1(gm1), 
                 err_rc_gm1, 0., "rc from g-1");

    hope.isclose(s1.moment_inertia_from_center_gm1(gm1), 
                 s2.moment_inertia_from_center_gm1(gm1), 
                 err_mi_gm1, 0., "mi from g-1");

    hope.isclose(s1.lambda_tidal_from_center_gm1(gm1), 
                 s2.lambda_tidal_from_center_gm1(gm1), 
                 err_lt_gm1, 0., "lt from g-1");
  }

  auto img = linear_spacing(rg_mg.min(), rg_mg.max(), nsamp);
  for (real_t mg : img) 
  {
    assert(s1.contains_grav_mass(mg));
    assert(s2.contains_grav_mass(mg));
    
    if (!hope.isclose(s1.center_gm1_from_grav_mass(mg), 
                 s2.center_gm1_from_grav_mass(mg), 
                 err_gm1_mg, 0., "g-1 from mg"))
    {
      hope.postmortem(str(boost::format(
        "mg = %.15e = (1 + %.15e) mg_max ") 
        % mg % (mg/rg_mg.max() - 1)));
    }
                 
    hope.isclose(s1.bary_mass_from_grav_mass(mg), 
                 s2.bary_mass_from_grav_mass(mg), 
                 err_m_mg, 0., "mb from mg");

    hope.isclose(s1.circ_radius_from_grav_mass(mg), 
                 s2.circ_radius_from_grav_mass(mg), 
                 err_rc_mg, 0., "rc from mg");

    hope.isclose(s1.moment_inertia_from_grav_mass(mg), 
                 s2.moment_inertia_from_grav_mass(mg), 
                 err_mi_mg, 0., "mi from mg");

    hope.isclose(s1.lambda_tidal_from_grav_mass(mg), 
                 s2.lambda_tidal_from_grav_mass(mg), 
                 err_lt_mg, 0., "lt from mg");
  }
  
  return hope;
}

auto get_eos_by_name(std::string s, units u)
->eos_barotr
{
    std::string eos_path{ std::string(PATH_TOV_EOS) + "/" 
                           + s + ".eos.h5" };
    return load_eos_barotr(eos_path, u);
}

bool test_branch_file_io(const star_branch seq, 
                         std::size_t nsamp=500, real_t tol=1e-13)
{
  char tmpn[L_tmpnam];
  {
    auto gotf{ std::tmpnam(tmpn) }; 
    assert(gotf);
  }
  
  save_star_branch(tmpn, seq);
  auto seq2 = load_star_branch(tmpn);
  std::remove(tmpn);
  
  //~ const real_t tol{ 1e-13 };
  failcount hope("Loading and saving star sequence branch works");

  hope(compare_star_seqs(seq, seq2, nsamp, tol, tol, tol, tol, 
                         tol, tol, tol, tol, tol, tol),
       "Star sequence still same after recovery from file");
  
  return hope;     
}

BOOST_AUTO_TEST_CASE( test_tovseq_acc )
{
  failcount hope("Computing TOV sequences for various EOS yields "
                 "same result as the reference data.");

  std::string eoslist[] = {
    "H4_Read_PP", 
    "H4_Read_PP.spline", 
    "WFF1_Read_PP", 
    "APR4_Read_PP", 
    "MPA1_Read_PP",
    "MS1_Read_PP"
    };
    
  auto u{ units::geom_solar() };
  
  for (auto s : eoslist) 
  {
    
    auto eos{ get_eos_by_name(s, u) };
    const tov_acc_simple acc{1e-8, 1e-6};
    auto seq = make_tov_branch_stable(eos, acc);

    std::string ref_path{ std::string(PATH_TOV_REF) + "/ref_tovseq_" 
                           + s + ".h5" };

    auto refseq = load_star_branch(ref_path, u);
    
    std::size_t nsamp{500};
    
    real_t tol_gm1{1e-12};
    real_t tol_mg{5e-10}; // Lower accuracy at maximum mass
    real_t err_m_gm1{tol_gm1}, err_rc_gm1{tol_gm1}, err_mi_gm1{tol_gm1}; 
    real_t err_lt_gm1{tol_gm1}, err_gm1_mg{tol_mg};
    real_t err_m_mg{tol_mg}, err_rc_mg{tol_mg}, err_mi_mg{tol_mg}; 
    real_t err_lt_mg{tol_mg};
    real_t err_rg{1e-12};
           
    if (!hope(compare_star_seqs(seq, refseq, nsamp,
                err_m_gm1, err_rc_gm1, err_mi_gm1, err_lt_gm1, 
                err_gm1_mg, err_m_mg, err_rc_mg, err_mi_mg, err_lt_mg,
                err_rg), 
         "comparing seqs for one EOS"))
    {
      hope.postmortem(str(boost::format("EOS = %s") %  s));
    }
    //~ save_star_branch(seq_path, seq);
  }
  
}

BOOST_AUTO_TEST_CASE( test_tovseq_file )
{
  failcount hope("Saving and loading the reference TOV sequences "
                 "works correctly");

  std::string eoslist[] = {
    "H4_Read_PP", 
    "H4_Read_PP.spline", 
    "WFF1_Read_PP", 
    "APR4_Read_PP", 
    "MPA1_Read_PP",
    "MS1_Read_PP"
    };
    
  auto u{ units::geom_solar() };
  
  for (auto s : eoslist) 
  {
    std::string ref_path{ std::string(PATH_TOV_REF) + "/ref_tovseq_" 
                           + s + ".h5" };
    auto refseq = load_star_branch(ref_path, u);
    
    hope(test_branch_file_io(refseq), 
         str(boost::format("File IO for sequence = %s") %  s));
  }
}
    
    
    
    
