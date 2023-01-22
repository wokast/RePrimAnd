#define BOOST_TEST_MODULE EOS

#include <cstdio>
#include <boost/test/unit_test.hpp>
#include "test_utils.h"
#include <boost/format.hpp>

#include "test_config.h"

#include "unitconv.h"
#include "eos_thermal.h"
#include "eos_thermal_file.h"
#include "eos_idealgas.h"
#include "eos_barotropic.h"
#include "eos_barotr_file.h"
#include "eos_barotr_poly.h"
#include "eos_barotr_spline.h"
#include "eos_hybrid.h"
#include "interpol.h"

#include "eos_data_ms1.h"

using boost::format;

using namespace std;
using namespace EOS_Toolkit;


bool check_eos_barotr_pw(const eos_barotr& eos, real_t rho, real_t eps,
                         real_t p, real_t cs, real_t gm1, 
                         real_t err_rho, real_t err_rhoe, real_t err_p,
                         real_t err_cs, real_t err_gm1) 
{
  failcount hope("Evaluating pressure, soundspeed, internal energy "
                 "from either mass density or pseudo-enthalpy");                    
                 
  if (eos.range_gm1().contains(gm1)) {
    auto s1 = eos.at_gm1(gm1);
    if (hope(s1.valid(), "Can evaluate EOS at valid g-1")) {
    
      hope.isclose(rho, s1.rho(), err_rho, 0., "rho from gm1");
      hope.isclose(1. + eps, 1. + s1.eps(),  err_rhoe, 0, 
                   "eps+1 from gm1");
      hope.isclose(p, s1.press(), err_p, 0, "P from gm1");
      hope.isclose(cs, s1.csnd(), 0., err_cs, "csnd from gm1");
    }
  }
  if (eos.range_rho().contains(rho)) {
    auto s2 = eos.at_rho(rho);
    if (hope(s2.valid(), "Can evaluate EOS at valid rho")) {  
      hope.isclose(gm1, s2.gm1(), err_gm1, 0., "gm1 from rho");
      hope.isclose(1.+eps, 1.+s2.eps(), err_rhoe, 0, "eps+1 from rho");
      hope.isclose(p, s2.press(), err_p, 0, "P from rho");
      hope.isclose(cs, s2.csnd(), 0., err_cs,  "csnd from rho");      
    }
  }
  if (!hope) {
    hope.postmortem(str(format("rho = %.15e, g-1 = %.15e") 
                         % rho % gm1));
  }
  
  return hope;

}


bool compare_eos_barotr(eos_barotr eos1, eos_barotr eos2, 
                        real_t rho0, real_t rho1, std::size_t nsamp,
                        real_t err_rho, real_t err_rhoe, real_t err_p,
                        real_t err_cs, real_t err_gm1)  
{
  failcount hope("Check that two EOS objects agree");
  
  auto irho = log_spacing(rho0, rho1, nsamp);
  for (real_t rho : irho) 
  {
    auto s1 = eos1.at_rho(rho);
    auto s2 = eos2.at_rho(rho);
    assert(s1);
    assert(s2);
    hope.isclose(s1.gm1(), s2.gm1(), err_gm1, 0., "gm1 from rho");
    hope.isclose(s1.press(), s2.press(), err_p, 0., "p from rho");
    hope.isclose(s1.csnd(), s2.csnd(), err_cs, 0., "csnd from rho");
    hope.isclose(1.+s1.eps(), 1.+s2.eps(), err_rhoe, 0., 
                 "1+eps from rho");
    
    real_t gm1 = s1.gm1();

    auto t1 = eos1.at_gm1(gm1);
    auto t2 = eos2.at_gm1(gm1);
    assert(t1);
    assert(t2);
    hope.isclose(t1.rho(), t2.rho(), err_rho, 0., "rho from gm1");
    hope.isclose(t1.press(), t2.press(), err_p, 0., "p from gm1");
    hope.isclose(t1.csnd(), t2.csnd(), err_cs, 0., "csnd from gm1");
    hope.isclose(1.+t1.eps(), 1.+t2.eps(), err_rhoe, 0., 
                 "1+eps from gm1");
    
  }
  return hope;
}

BOOST_AUTO_TEST_CASE( test_eos_lookup_table )
{
  failcount hope("Lookup tables used in barotropic EOS "
                 "works");
                 
  const std::size_t npts = 400;           
  const real_t eqtol = 1e-14;
      
  auto f1 = [] (real_t x) {
    return x*(x + 10.) + sin(x); 
  };
  // |second derivative| <= 3
  
  lookup_table t1{f1, {-5., 100.}, npts};
  
  auto ix = linear_spacing(t1.range_x().min(),
                           t1.range_x().max(), npts-1);
  for (real_t x : ix) {
    hope.isclose(f1(x), t1(x), eqtol, 0,
                 "Linear lookup at tabulated points");
  }

  // expect error ~ f'' h^2
  const real_t tol   = pow(t1.range_x().length() / npts, 2);

  auto ix1 = linear_spacing(t1.range_x().min(),
                           t1.range_x().max(), 10*npts);
  for (real_t x : ix1) {
    hope.isclose(f1(x), t1(x), 0, tol,
                 "Linear lookup between tabulated points");
  }

  auto f2 = [] (real_t x) {
    return 1.23456 * x + 789; 
  };
  
  lookup_table t2{f2, {-5., 100.}, npts};
  
  auto ix2 = linear_spacing(t2.range_x().min(),
                            t2.range_x().max(), 10 * npts);
  for (real_t x : ix2) {
    hope.isclose(f2(x), t2(x), eqtol, 0,
                 "Linear lookup of linear function");
  }

  const int mags = 10;
  const real_t mags10 = pow(10., -mags);
  
  auto f3 = [] (real_t x) {
    return 2.4*pow(log(x),2) + 3.0; 
  };

  lookup_table_magx m3{f3, {1e5*mags10, 1e5}, npts, mags};
  
  auto ix3 = log_spacing(m3.range_x().min(),
                         m3.range_x().max(), npts-1);
  for (real_t x : ix3) {
    hope.isclose(f3(x), m3(x), eqtol, 0,
                 "Log lookup at sample points");
  }
  
  lookup_table::range_t rg4{1e-7, 1e5};
  const real_t ofs = (rg4.max()*mags10 - rg4.min()) / (1.0-mags10);
  
  auto f4 = [ofs] (real_t x) {
    return 2*log(x+ofs) + 3.0; 
  };

  lookup_table_magx m4{f4, rg4, npts, mags};
  
  auto ix4 = log_spacing(rg4.min()+ofs, rg4.max()+ofs, 10 * npts);
  for (real_t x : ix4) {
    real_t z = x - ofs;
    hope.isclose(f4(z), m4(z), 0, 10*eqtol,
                 "Log lookup of log function");
  }
  

}



BOOST_AUTO_TEST_CASE( test_eos_barotr_fail )
{
  failcount hope("barotropic EOS throws as expected when "
                 "called uninitialized or on invalid matter state");

  eos_barotr eos{};
  hope.dothrow("invalid_eos.is_isentropic()", 
               [&] () {eos.is_isentropic();});
  hope.dothrow("invalid_eos.is_zero_temp()", 
               [&] () {eos.is_zero_temp();});
  hope.dothrow("invalid_eos.has_temp()", 
               [&] () {eos.has_temp();});
  hope.dothrow("invalid_eos.has_efrac()", 
               [&] () {eos.has_efrac();});
  hope.dothrow("invalid_eos.range_rho()", 
               [&] () {eos.range_rho();});
  hope.dothrow("invalid_eos.range_gm1()", 
               [&] () {eos.range_gm1();});
  hope.dothrow("invalid_eos.minimal_h()", 
               [&] () {eos.minimal_h();});
  hope.dothrow("invalid_eos.is_rho_valid()", 
               [&] () {eos.is_rho_valid(0.0);});
  hope.dothrow("invalid_eos.is_gm1_valid()", 
               [&] () {eos.is_gm1_valid(0.0);});
  hope.dothrow("invalid_eos.at_rho()", 
               [&] () {eos.at_rho(0.0);});
  hope.dothrow("invalid_eos.at_gm1()", 
               [&] () {eos.at_gm1(0.0);});
  
  const real_t rhomax = 10.;
  eos = make_eos_barotr_poly(2., 0.1, rhomax);

  hope.nothrow("valid_eos.range_rho()", 
               [&] () {eos.range_rho();});
  hope.nothrow("valid_eos.at_rho(valid_rho)", 
               [&] () {eos.at_rho(rhomax/2);});
  hope.nothrow("valid_eos.at_rho(invalid_rho)", 
               [&] () {eos.at_rho(rhomax*2);});

  auto s = eos.at_rho(1.001*rhomax);

  hope(!s, "bool(state(invalid_rho)) == false");

  hope.dothrow("invalid_state.gm1()",   [&] () {s.gm1();}   );
  hope.dothrow("invalid_state.rho()",   [&] () {s.rho();}   );
  hope.dothrow("invalid_state.press()", [&] () {s.press();} );
  hope.dothrow("invalid_state.eps()",   [&] () {s.eps();}   );
  hope.dothrow("invalid_state.hm1()",   [&] () {s.hm1();}   );
  hope.dothrow("invalid_state.csnd()",  [&] () {s.csnd();}  );
  hope.dothrow("invalid_state.temp()",  [&] () {s.temp();}  );
  hope.dothrow("invalid_state.ye()",    [&] () {s.ye();}    );

  auto s2 = eos.at_rho(0.999*rhomax);

  hope(s2.valid(), "bool(state(valid_rho)) == true");

  hope.nothrow("valid_state.gm1()",   [&] () {s2.gm1();}   );
  hope.nothrow("valid_state.rho()",   [&] () {s2.rho();}   );
  hope.nothrow("valid_state.press()", [&] () {s2.press();} );
  hope.nothrow("valid_state.eps()",   [&] () {s2.eps();}   );
  hope.nothrow("valid_state.hm1()",   [&] () {s2.hm1();}   );
  hope.nothrow("valid_state.csnd()",  [&] () {s2.csnd();}  );
  hope.nothrow("valid_state.temp()",  [&] () {s2.temp();}  );

  hope.dothrow("EOS without Ye: valid_state.ye()",    
               [&] () {s2.ye();}    );

  auto s3 = eos.at_gm1(0.999*eos.range_gm1().max());

  hope(s3.valid(), "bool(state(valid_gm1)) == true");
  hope.nothrow("state(valid_gm1)).press()", [&] () {s3.press();} );
  
  auto s4 = eos.at_gm1(1.001*eos.range_gm1().max());

  hope(!s4.valid(), "bool(state(invalid_gm1)) == false");
  hope.dothrow("state(invalid_gm1)).press()", [&] () {s4.press();} );

}

BOOST_AUTO_TEST_CASE( test_eos_barotr_table )
{
  failcount hope("Tabulated barotropic EOS loaded from file "
                 "evaluates to correct values over full range");
                 
  auto u = units::geom_solar();
  eos_barotr eos = load_eos_barotr(PATH_EOS, u);
  
  
  for(unsigned int k=1; k<eos_data_ms1.size()-1; ++k) {
    const auto& m = eos_data_ms1[k]; 
    const auto& l = eos_data_ms1[k-1]; 
    const auto& r = eos_data_ms1[k+1];

    real_t rho = m[0] / u.density();
    real_t eps = m[1];
    real_t p   = m[2] / u.pressure();
    real_t cs  = sqrt(m[3]);
    real_t gm1 = m[4];
    
    const bool jump{ 
      fabs((r[3] - l[3]) * m[4] / (m[3] * (r[4]-l[4]))) > 10
    };

    real_t err_cs  = jump ? 5e-3 : 4e-4;    //EOS has jumps in csnd
    real_t err_gm1 = 5e-5;  
    real_t err_rho = 3e-4;
    real_t err_rhoe = 1e-5;
    real_t err_p   = 1e-4;
  
    
    hope(check_eos_barotr_pw(eos, rho, eps, p, cs, gm1, 
                             err_rho, err_rhoe, err_p, err_cs, err_gm1),
         "Evaluate EOS at single point");
    
  }
  
                 
}


BOOST_AUTO_TEST_CASE( test_eos_barotr_pwpoly_file )
{
  failcount hope("Loading EOS from file not failing");
                 
  auto u = units::geom_solar();
  
  hope.nothrow("loading PW poly",   
    [&] () {
      auto eos = load_eos_barotr(PATH_EOS_PP, u);
    } 
  );

  hope.nothrow("loading tabulated",   
    [&] () {
      auto eos = load_eos_barotr(PATH_EOS, u);
    } 
  );

}

std::string get_temp_filename()
{
  char tmpn[L_tmpnam];
  {
    auto gotf = std::tmpnam(tmpn); 
    assert(gotf);
  }
  return tmpn;  
}

  
bool test_eos_file_io(const eos_barotr eos, real_t rho0, 
                      std::size_t nsamp)
{
  auto tmpn = get_temp_filename();
  save_eos_barotr(tmpn, eos);
  auto eos2 = load_eos_barotr(tmpn);
  std::remove(tmpn.c_str());
  
  real_t tol = 1e-13;
  failcount hope("Loading and saving EOS works");

  hope(compare_eos_barotr(eos, eos2, 
                     rho0, eos.range_rho().max()*0.99, nsamp,
                     tol, tol, tol, tol, tol),
       "EOS still same after recovery from file");
 
  
  hope.isclose(eos.range_gm1().max(), eos2.range_gm1().max(), 
               tol, 0, "max valid range g-1");
  hope.isclose(eos.range_rho().max(), eos2.range_rho().max(), 
               tol, 0, "max valid range rho");
      
  return hope;
}
  

BOOST_AUTO_TEST_CASE( test_eos_spline )
{
  failcount hope("Spline EOS accurate");
                 
  auto u = units::geom_solar();
 
  const real_t n_poly      = 1;
  const real_t rmd_poly    = 6.176e+18 / u.density();
  const real_t rhomax_poly = 1E19 / u.density();
  auto eos_poly = make_eos_barotr_poly(n_poly, rmd_poly, rhomax_poly);
 
  const std::size_t pts_per_mag = 400;
  const int mags = 10;
  const real_t rho_max = 1e19 / u.density();
  const eos_barotr::range rg_spl{rho_max / pow(10,mags), rho_max};
  auto eos_spline = make_eos_barotr_spline(eos_poly, rg_spl,  
                                           n_poly, pts_per_mag);
                                           
  real_t err_cs  = 2e-6;
  real_t err_gm1 = 1e-6;  
  real_t err_rho = 1e-6;
  real_t err_rhoe = 2e-5;
  real_t err_p   = 1e-6;
  
  std::size_t nsamp = 3*(mags*pts_per_mag+2);
  
  hope(compare_eos_barotr(eos_poly, eos_spline, 
                     rg_spl.min()/100, rg_spl.max()*0.99, nsamp,
                     err_rho, err_rhoe, err_p, err_cs, err_gm1),
       "Compare spline representation of polytropic EOS to original");
  
  
  hope(test_eos_file_io(eos_poly, rg_spl.min()/100, nsamp),
       "Saving+loading polytropic EOS works");
  hope(test_eos_file_io(eos_spline, rg_spl.min()/100, nsamp),
       "Saving+loading spline EOS works");
}
  

BOOST_AUTO_TEST_CASE( test_eos_spline_fromtab )
{
  failcount hope("Spline EOS from table accurate");
                 
  auto u = units::geom_solar();
  std::vector<real_t> vgm1, vrho, veps, vpress, vcsnd, 
                            vtemp, vefrac;
  
  for(unsigned int k=1; k<eos_data_ms1.size()-1; ++k) {
    const auto& m = eos_data_ms1[k]; 

    vrho.push_back(m[0] / u.density());
    veps.push_back(m[1]);
    vpress.push_back(m[2] / u.pressure());
    vcsnd.push_back(sqrt(m[3]));
    vgm1.push_back(m[4]);
  };
  const real_t n_poly = 1.7115960633290546;
  const std::size_t pts_per_mag = 500;
  const int mags = 10;
  
  
  real_t rho0 = vrho.back() / pow(10, mags);
  std::size_t ij = 0;
  for (; ij<vrho.size(); ++ij) if (vrho[ij] > rho0) break;
  rho0 = vrho[ij];

  auto eos = make_eos_barotr_spline(
    vgm1, vrho, veps, vpress, vcsnd, vtemp, vefrac, true,            
    {rho0, vrho[vrho.size()-2]}, n_poly, u, pts_per_mag);


  const real_t gcorr{ 
    (vgm1[ij] - eos.gm1_at_rho(vrho[ij])) / (1.0+vgm1[ij])
  };


  for (auto m : eos_data_ms1) {
    real_t rho = m[0] / u.density();
    
    if (rho <= rho0) continue;
    
    real_t eps = m[1];
    real_t p   = m[2] / u.pressure();
    real_t cs  = sqrt(m[3]);
    real_t gm1raw = m[4];
    
    real_t gm1 = gm1raw - gcorr * (1.0 + gm1raw);
    


    real_t err_cs  = 1e-2;    //EOS has jumps in csnd
    real_t err_gm1 = 1e-5;  
    real_t err_rho = 3e-3;    //EOS has steep gradients and kinks
    real_t err_rhoe = 1e-5;
    real_t err_p   = 2e-4;
   
    hope(check_eos_barotr_pw(eos, rho, eps, p, cs, gm1, 
                             err_rho, err_rhoe, err_p, err_cs, err_gm1),
         "Evaluate EOS at single point");
    
  }
  
  auto eos2 = make_eos_barotr_spline(
                  vrho, veps, vpress, vcsnd, vtemp, vefrac, true,            
                  {rho0, vrho[vrho.size()-2]}, n_poly, u, pts_per_mag);

  hope(compare_eos_barotr(eos, eos2, 
                          rho0, eos2.range_rho().max()*0.99, 400,
                          6e-3,  //err_rho, 
                          5e-4,  //err_rhoe, 
                          1e-3,  //err_p,
                          3e-2,  //err_cs, 
                          5e-4), //err_gm1) ,
       "Creating EOS from samples with and without gm1 equivalent");

  
}




BOOST_AUTO_TEST_CASE( test_eos_thermal_save )
{
  failcount hope("Can load and save thermal EOS");
  
  auto u = units::geom_solar(); 
  real_t eps_max{ 1e2 };
  real_t rho_max{ 0.1 };
  
  auto eos1 = make_eos_idealgas(1.0, eps_max, rho_max);
  
  auto fn1 = get_temp_filename();
  hope.nothrow("Can save ideal gas EOS", [&] () {
    save_eos_thermal(fn1, eos1);
  });
  hope.nothrow("Can load ideal gas EOS", [&] () {
    auto eos2 = load_eos_thermal(fn1, u);
  });
  std::remove(fn1.c_str());
  
  auto eos3 = load_eos_barotr(PATH_EOS_PP, u);
  auto eos4 = make_eos_hybrid(eos3, 1.8, eps_max, rho_max);

  auto fn2 = get_temp_filename();
  hope.nothrow("Can save hybrid (based on pw. poly) EOS", [&] () {
    save_eos_thermal(fn2, eos4);
  });
  hope.nothrow("Can load hybrid (based on pw. poly) EOS", [&] () {
    auto eos5 = load_eos_thermal(fn2, u);
  });
  std::remove(fn2.c_str());
  
}




