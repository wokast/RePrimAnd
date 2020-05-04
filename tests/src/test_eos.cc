#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE EOS

#include <boost/test/unit_test.hpp>
#include "test_utils.h"
#include <boost/format.hpp>

#include "test_config.h"

#include "unitconv.h"
#include "eos_thermal.h"
#include "eos_idealgas.h"
#include "eos_barotropic.h"
#include "eos_barotr_file.h"
#include "eos_barotr_poly.h"
#include "eos_hybrid.h"
#include "interpol.h"

#include "eos_data_ms1.h"

using boost::format;

using namespace std;
using namespace EOS_Toolkit;


bool check_eos_barotr_pw(const eos_barotr& eos, real_t rho, real_t eps,
                         real_t p, real_t cs, real_t gm1, 
                         real_t err_rho, real_t err_eps, real_t err_p,
                         real_t err_cs, real_t err_gm1) 
{
  failcount hope("Evaluating pressure, soundspeed, internal energy "
                 "from either mass density or pseudo-enthalpy");                    
                 
  {
    auto s1 = eos.at_gm1(gm1);
    if (hope(s1.valid(), "Can evaluate EOS at valid g-1")) {
    
      hope.isclose(rho, s1.rho(), err_rho, 0, "rho from gm1");
      hope.isclose(eps, s1.eps(), err_eps, 0, "eps from gm1");
      hope.isclose(p, s1.press(), err_p, 0, "P from gm1");
      hope.isclose(cs, s1.csnd(), err_cs, 0, "csnd from gm1");
    }
  }
  {
    auto s2 = eos.at_rho(rho);
    if (hope(s2.valid(), "Can evaluate EOS at valid rho")) {  
      hope.isclose(gm1, s2.gm1(), err_gm1, 0, "gm1 from rho");
      hope.isclose(eps, s2.eps(), err_eps, 0, "eps from rho");
      hope.isclose(p, s2.press(), err_p, 0, "P from rho");
      hope.isclose(cs, s2.csnd(), err_cs, 0, "csnd from rho");      
    }
  }
  if (!hope) {
    hope.postmortem(str(format("rho = %.15e, g-1 = %.15e") 
                         % rho % gm1));
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


BOOST_AUTO_TEST_CASE( test_eos_splines )
{
  failcount hope("Spline interpolation used to load tabulated EOS "
                 "from file works");
                 
  const std::size_t npts = 400;           
  const real_t eqtol = 1e-14;
      
  auto f1 = [] (real_t x) {
    return x*(x + 10.) + sin(x); 
  };
  // |second derivative| <= 3
  
  std::vector<real_t> x1,y1,x2,y2;
  
  
  auto ix = linear_spacing(-8., 123.4, npts-1);

  for (real_t x : ix) {
    x1.push_back(x);
    y1.push_back(f1(x)); 
  }
  
  cspline_mono s1{x1,y1};
  
  for (real_t x : ix) {
    hope.isclose(f1(x), s1(x), eqtol, 0,
                 "Spline approximation at sample points");
  }

  // expect error ~ f''' h^3
  const real_t tol   = pow(s1.range_x().length() / npts, 3);

  auto ix1 = linear_spacing(s1.range_x().min(),
                            s1.range_x().max(), 10*npts);
  for (real_t x : ix1) {
    hope.isclose(f1(x), s1(x), 0, tol,
                 "Spline approximation between sample points");
  }

  auto f2 = [] (real_t x) {
    real_t z = x + 10;
    return 1.23456 * z*z + 789; 
  };
  
  for (real_t x : ix) {
    x2.push_back(x);
    y2.push_back(f2(x)); 
  }

  cspline_mono s2{x2,y2};
  
  const real_t margin = s2.range_x().length() / (npts-1);
  auto ix2 = linear_spacing(s2.range_x().min() + margin,
                            s2.range_x().max() - margin, 10 * npts);
  for (real_t x : ix2) {
    hope.isclose(f2(x), s2(x), eqtol, 0,
                 "spline approximation for quadratic function");
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

    real_t err_cs  = jump ? 5e-2 : 4e-4;    //EOS has jumps in csnd
    real_t err_gm1 = 5e-5;  
    real_t err_rho = 3e-4;
    real_t err_eps = 5e-5;
    real_t err_p   = 1e-4;
  
    
    hope(check_eos_barotr_pw(eos, rho, eps, p, cs, gm1, 
                             err_rho, err_eps, err_p, err_cs, err_gm1),
         "Evaluate EOS at single point");
    
  }
  
                 
}
