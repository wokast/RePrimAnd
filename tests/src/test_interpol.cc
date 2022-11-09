#define BOOST_TEST_MODULE INTERPOL

#include <boost/test/unit_test.hpp>
#include "test_utils.h"
#include <boost/format.hpp>

#include "test_config.h"

#include "unitconv.h"
#include "interpol.h"
#include "interpol_regspl.h"
#include "interpol_logspl.h"
#include "interpol_pchip_spline.h"
#include "interpol_linear.h"



using namespace std;
using namespace EOS_Toolkit;


                         

BOOST_AUTO_TEST_CASE( test_interp_linear )
{
  failcount hope("Test linear interpolation methods");
                 
  const std::size_t npts = 400;           
  const real_t eqtol = 1e-14;
      
  auto f1 = [] (real_t x) {
    return x*(x + 10.) + sin(x); 
  };
  
  auto t1 = make_interpol_reglin(f1, {-5., 100.}, npts);
  
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
  auto t2 = make_interpol_reglin(f2, {-5., 100.}, npts);

  
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

  auto m3 = make_interpol_loglin(f3, {1e5*mags10, 1e5}, npts);
  auto ix3 = log_spacing(m3.range_x().min(),
                         m3.range_x().max(), npts-1);
  for (real_t x : ix3) {
    hope.isclose(f3(x), m3(x), eqtol, 0,
                 "Log lookup at sample points");
  }
  
  interpolator::range_t rg4{1e-7, 1e5};
  
  auto f4 = [] (real_t x) {
    return 2*log(x) + 3.0; 
  };

  auto m4 = make_interpol_loglin(f4, rg4, npts);
  
  auto ix4 = log_spacing(rg4.min(), rg4.max(), 10 * npts);
  for (real_t x : ix4) {
    hope.isclose(f4(x), m4(x), eqtol, eqtol,
                 "Log lookup of log function");
  }
  
}


BOOST_AUTO_TEST_CASE( test_interp_reg_spline )
{
  failcount hope("Regular monotonic spline works");
       
  
  const std::size_t npts = 400;           
  const real_t eqtol = 1e-14;
  const real_t eqtol_abs = 1e-12;
  
  auto f1 = [] (real_t x) {
    return x*(x + 10.) + sin(x); 
  };
  // |second derivative| <= 3
  
  auto t1 = make_interpol_regspl(f1, {-5., 100.}, npts);
  
  auto ix = linear_spacing(t1.range_x().min(),
                           t1.range_x().max(), npts-1);
  for (real_t x : ix) {
    hope.isclose(f1(x), t1(x), eqtol, eqtol_abs,
                 "Regular spline at tabulated points");
  }

  const real_t tol   = 10*pow(t1.range_x().length() / npts, 2);

  auto ix1 = linear_spacing(t1.range_x().min(),
                           t1.range_x().max(), 10*npts);
  for (real_t x : ix1) {
    hope.isclose(f1(x), t1(x), 0, tol,
                 "Regular spline between tabulated points");
  }

  auto f2 = [] (real_t x) {
    return 1.23456 * x + 789; 
  };
  
  auto t2 = make_interpol_regspl(f2, {-5., 100.}, npts);
  
  auto ix2 = linear_spacing(t2.range_x().min(),
                            t2.range_x().max(), 10 * npts);
  for (real_t x : ix2) {
    hope.isclose(f2(x), t2(x), eqtol, eqtol_abs,
                 "Regular spline of linear function");
  }

  const int mags = 10;
  const real_t mags10 = pow(10., -mags);
  
  auto f3 = [] (real_t x) {
    return 2.4*pow(log(x),2) + 3.0; 
  };

  auto m3 = make_interpol_logspl(f3, {1e5*mags10, 1e5}, npts);
  
  auto ix3 = log_spacing(m3.range_x().min(),
                         m3.range_x().max(), npts-1);
  for (real_t x : ix3) {
    hope.isclose(f3(x), m3(x), eqtol, eqtol_abs,
                 "Log spline at sample points");
  }
  
  interpolator::range_t rg4{1e-7, 1e5};
  
  auto f4 = [] (real_t x) {
    return 2*log(x) + 3.0; 
  };
  
  auto m4 = make_interpol_logspl(f4, rg4, npts);
  
  auto ix4 = log_spacing(rg4.min(), rg4.max(), 10 * npts);
  for (real_t x : ix4) {
    hope.isclose(f4(x), m4(x), eqtol, eqtol_abs,
                 "Log spline of log function");
  }
  

  interpolator::range_t rg5{1e-10, 1e10};
  
  auto f5 = [] (real_t x) {
    return x*x; 
  };

  auto m5 = make_interpol_llogspl(f5, rg5, npts);
  
  auto ix5 = log_spacing(rg5.min(), rg5.max(), 10 * npts);
  for (real_t x : ix5) {
    hope.isclose(f5(x), m5(x), 10*eqtol, 0,
                 "Log-log spline interpolation high dyn range");
  }
}


BOOST_AUTO_TEST_CASE( test_interp_pchip_spline )
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
  
  auto s1 = make_interpol_pchip_spline(x1, y1);

  
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

  auto s2 = make_interpol_pchip_spline(x2, y2);
  
  const real_t margin = s2.range_x().length() / (npts-1);
  auto ix2 = linear_spacing(s2.range_x().min() + margin,
                            s2.range_x().max() - margin, 10 * npts);
  for (real_t x : ix2) {
    hope.isclose(f2(x), s2(x), eqtol, 0,
                 "spline approximation for quadratic function");
  }

}


  
  
  
  
