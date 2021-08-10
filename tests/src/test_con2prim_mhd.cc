#define BOOST_TEST_MODULE con2prim
#include <boost/test/unit_test.hpp>
#include <limits>

#include "test_config.h"
#include "unitconv.h"
#include "test_con2prim_mhd.h"

#include "eos_thermal.h"
#include "eos_idealgas.h"
#include "eos_thermal_file.h"
#include "eos_hybrid.h"

using boost::format;
using boost::str;

using namespace std;
using namespace EOS_Toolkit;


real_t w_from_z(real_t z) {
  return sqrt(1.0+z*z);
}

real_t v_from_z(real_t z) {
  return z/w_from_z(z);
}


void test_con2prim_mhd::setup_prim_cons(
        prim_vars_mhd& pv, cons_vars_mhd& cv, 
        real_t rho, real_t eps, real_t ye, 
        real_t z, real_t b, int vdim, int bdim) const
{
  real_t wl0 = sqrt(1.0 + z*z);
  real_t v   = z / wl0;
  
  //~ eos_thermal::status errs;
  real_t press0 = 0.0;
  if (rho > 0) {
    auto s = eos.at_rho_eps_ye(rho, eps, ye);
    press0 = s.press();    
  }
  sm_vec3u vel0 = ZERO;
  vel0(vdim)    = v;
  //~ real_t wl0    = 1.0 / sqrt(1.0 - v*v);
  sm_vec3u B0   = ZERO;
  B0(bdim)      = b*sqrt(rho*wl0);
  sm_vec3l El   = g.cross_product(B0, vel0);
  sm_vec3u E0   = g.raise(El);

  pv            = prim_vars_mhd(rho, eps, ye, press0, 
                                vel0, wl0, E0, B0);
  cv.from_prim(pv, g);
}


bool test_con2prim_mhd::check_c2p(prim_vars_mhd& pv, 
                 cons_vars_mhd& cv, 
                 con2prim_mhd::report& rep) const
{
  failcount hope("Calling C2P works");
  
  hope.nothrow("C2P", cv2pv, pv, cv, g, rep);
  if (!hope) return hope;

  if (rep.failed()) {
    hope(check_isnan(cv,pv), 
         "C2P failure implies all results set to NAN");
    hope(rep.adjust_cons, 
         "C2P failure implies report adjust conserved flag set");
    hope(!rep.set_atmo, 
         "C2P report atmo set flag implies not failed");      
  } else {
    hope(check_isfinite(cv, pv), 
         "All C2P results finite");
    hope(check_vel_w(pv), 
         "C2P computes consistent velocity vector and Lorentz factor");
    hope(eos.is_rho_eps_ye_valid(pv.rho, pv.eps, pv.ye),
         "C2P results inside EOS validity range");
  }
  if (rep.set_atmo) {
    hope(rep.adjust_cons, 
         "C2P report set atmo flag implies adjust conserved flag set");
  }
  
  return hope;
}


bool test_con2prim_mhd::check_same(const prim_vars_mhd&pv0, 
                const prim_vars_mhd& pv1) const
{
  failcount hope("Primitive variables agree up to rounding errors");

  const real_t tol = 1e-13;
  
  hope.isclose(pv1.rho, pv0.rho, tol, 0, "rho");
  hope.isclose(pv1.eps, pv0.eps, tol, 0, "eps");
  hope.isclose(pv1.ye, pv0.ye, tol, 0, "ye");
  hope.isclose(pv1.press, pv0.press, tol, 0, "press");
  hope.isclose(pv1.w_lor, pv0.w_lor, tol, 0, "w_lor");
  hope(check_isclose(pv1.vel, pv0.vel, tol, 0), "E");
  hope(check_isclose(pv1.E, pv0.E, tol, 0), "E");
  hope(check_isclose(pv1.B, pv0.B, tol, 0), "B");

  return hope;
}

bool test_con2prim_mhd::check_same(const cons_vars_mhd& cv0, 
                const cons_vars_mhd& cv1) const
{
  failcount hope("Conserved variables agree up to rounding errors");

  const real_t tol{ 1e-13 };
    
  hope.isclose(cv1.dens, cv0.dens, tol, 0, "dens");
  hope.isclose(cv1.tau, cv0.tau, tol, 0, "tau");
  hope.isclose(cv1.tracer_ye, cv0.tracer_ye, tol, 0, "tracer_ye");
  hope(check_isclose(cv1.scon, cv0.scon, tol, 0), "scon");
  hope(check_isclose(cv1.bcons, cv0.bcons, tol, 0), "bcons");
  
  return hope;
} 


bool test_con2prim_mhd::compare_prims(
        const prim_vars_mhd& pv0, 
        const prim_vars_mhd& pv1) const
{
  failcount hope("Primitive variables agree within "
                 "expected C2P tolerances");
    
  
  hope((pv1.rho > 0) && (pv0.rho > 0), "Densities positive");
  hope((pv1.eps > -1) && (pv0.eps > -1), "Energy densities positive");
  

  const real_t vsqr0 = g.norm2(pv0.vel);
  const real_t vsqr1 = g.norm2(pv1.vel);
  const real_t vsqr  = max(vsqr0,vsqr1);
  const real_t v  = sqrt(vsqr);
  const real_t dv = sqrt(g.norm2(pv1.vel - pv0.vel));
  
  const real_t w  = max(pv0.w_lor, pv1.w_lor);
  const real_t wsqr = w*w;
  
  const real_t eps  = (pv0.eps + pv1.eps) / 2;
  const real_t deps = fabs(pv1.eps - pv0.eps);
  
  const real_t Esqr0 = g.norm2(pv0.E);
  const real_t Esqr1 = g.norm2(pv1.E);
  const real_t E  = sqrt(max(Esqr0,Esqr1));
  const real_t dE = sqrt(g.norm2(pv1.E - pv0.E));
  
  const real_t Bsqr0 = g.norm2(pv0.B);
  const real_t Bsqr1 = g.norm2(pv1.B);
  const real_t Bsqr  = std::max(Bsqr0, Bsqr1);
  const real_t B  = sqrt(Bsqr);
  const real_t dB = sqrt(g.norm2(pv1.B - pv0.B));

  const real_t a0 = pv0.press / (pv0.rho*(1.0 + pv0.eps));
  const real_t a1 = pv1.press / (pv1.rho*(1.0 + pv1.eps));
  const real_t a  = std::max(a0,a1);

  const real_t cs0  = eos.at_rho_eps_ye(pv0.rho, pv0.eps, pv0.ye).csnd();
  const real_t cs1  = eos.at_rho_eps_ye(pv1.rho, pv1.eps, pv1.ye).csnd();
  const real_t csqr = pow(std::max(cs0,cs1),2);
  
  const real_t dfmin = 1.0 - vsqr*csqr;


  const real_t acc_fudge = 10;
  const real_t acc_max { 
    acc_fudge * std::numeric_limits<real_t>::epsilon() 
  }; 
  
  const real_t acc{ cv2pv.get_acc() + acc_max * wsqr / dfmin }; 

  
  const real_t acc_eps {
    acc_max * std::max({w * (Bsqr / (pv0.rho * w)) / eps, 
                        vsqr*wsqr / eps, 
                        1.0})
  }; 
  const real_t acc_w { acc_max * wsqr };

  
  if ((pv0.press == 0)  || (pv1.press ==0)) {
    //special case, e.g. ideal Gas at T=0
    hope(pv0.press == pv1.press, "press equal");
  } 
  else {        
    const real_t f0 = cs0*cs0 * (1.0 + a0) / a0;
    const real_t f1 = cs1*cs1 * (1.0 + a1) / a1;
    const real_t f  = max(f0,f1);
    
    const real_t F = (pv0.eps/pv0.press) 
                       * eos.at_rho_eps_ye(pv0.rho, 
                                      pv0.eps, pv0.ye).dpress_deps();
    
    hope.isclose(pv0.press, pv1.press, f*acc + F*acc_eps, 0, 
                 "press");
  }  
  hope(dv <= v * (acc/wsqr + acc_max), 
       "v within tolerance");
       
  hope.isclose(pv1.w_lor, pv0.w_lor, vsqr * acc + acc_w, 0, "W");
  hope.isclose(pv1.rho, pv0.rho, vsqr * acc + acc_w, 0, "rho");
  
  hope(deps <= (1.0 + eps) * (a * vsqr * acc + acc_eps), 
       "specific total energy within tolerance");
  hope(dE <= E * (acc/wsqr + acc_max), 
       "E within tolerance");
  
  hope(dB <= B * acc_max, "B identical");
  hope.isclose(pv1.ye, pv0.ye, acc_max, 0, "ye");
  
  cons_vars_mhd cv0, cv1;
  cv0.from_prim(pv0, g);
  cv1.from_prim(pv1, g);
  compare_cons(cv0, cv1, vsqr);  
  
  return hope;
}

bool test_con2prim_mhd::compare_cons(
         const cons_vars_mhd& cv0, 
         const cons_vars_mhd& cv1, real_t vsqr) const
{
  failcount hope("Conserved variables agree within "
                 "expected C2P tolerances");
  
  const real_t acc{ cv2pv.get_acc() }; 
  const real_t rounding{ 1e-13 };
  
  hope.isclose(cv1.dens, cv0.dens, 0, rounding, "dens");
  hope.isclose(1. + cv1.tau, 1. + cv0.tau, 
               4 * vsqr * acc + rounding, 0, 
               "tau");
  hope.isclose(cv1.tracer_ye, cv0.tracer_ye, rounding, 0, 
               "tracer_ye");
  hope(check_isclose(cv1.scon, cv0.scon, acc, 0), "scon");
  hope(check_isclose(cv1.bcons, cv0.bcons, rounding, 0), "bcons");
  
  return hope;
}

bool test_con2prim_mhd::check_vel_w(const prim_vars_mhd& pv) const
{
  failcount hope("Primitive vars velocity vector and Lorentz factor"
                 "valid and consistent");
                 
  real_t reltol = 1e-14;
               
  const real_t vsqr = g.norm2(pv.vel);
  hope.isfinite(vsqr, "v^i v_i");  
  hope.isless(vsqr, 1.0, "v^i v_i");
  const real_t wsqr  = 1.0 / (1.0 - vsqr);
  const real_t w_lor = sqrt(wsqr);
  hope.isfinite(wsqr, "W^2(v^i v_i)");
  hope.isfinite(w_lor, "W(v^i v_i)");
  hope.isfinite(pv.w_lor, "prim.w_lor");

  hope.isclose(w_lor, pv.w_lor, wsqr*reltol, 0, "W");
         
  return hope;
}



bool test_con2prim_mhd::check_isnan(const cons_vars_mhd& cv) const
{
  failcount hope("All conserved vars set to NAN");
  
  hope.isnan(cv.dens, "dens"); 
  hope.isnan(cv.tau, "tau"); 
  hope.isnan(cv.tracer_ye, "tracer_ye");
  hope(check_isnan(cv.scon), "scon");  
  hope(check_isnan(cv.bcons), "bcons");
  
  return hope;
}

bool test_con2prim_mhd::check_isnan(const prim_vars_mhd& pv) const
{
  failcount hope("All primitive vars set to NAN");
  
  hope.isnan(pv.rho, "rho"); 
  hope.isnan(pv.eps, "eps"); 
  hope.isnan(pv.ye, "ye");
  hope.isnan(pv.press, "press");
  hope.isnan(pv.w_lor, "w_lor");
  hope(check_isnan(pv.vel), "vel is NAN"); 
  hope(check_isnan(pv.E), "E is NAN");
  hope(check_isnan(pv.B), "B is NAN");
  
  return hope;
}

bool test_con2prim_mhd::check_isnan(
          const cons_vars_mhd& cv,
          const prim_vars_mhd& pv) const
{
  return check_isnan(pv) && check_isnan(cv);
}



bool test_con2prim_mhd::check_isfinite(
          const cons_vars_mhd& cv) const
{
  failcount hope("All conserved vars finite");
  
  hope.isfinite(cv.dens, "dens"); 
  hope.isfinite(cv.tau, "tau"); 
  hope.isfinite(cv.tracer_ye, "tracer_ye");
  hope(check_isfinite(cv.scon), "scon finite");
  hope(check_isfinite(cv.bcons), "bcons finite");

  return hope;
}

bool test_con2prim_mhd::check_isfinite(
          const prim_vars_mhd& pv) const
{
  failcount hope("All primitive vars finite");
  
  hope.isfinite(pv.rho, "rho"); 
  hope.isfinite(pv.eps, "eps"); 
  hope.isfinite(pv.ye, "ye");
  hope.isfinite(pv.press, "press");
  hope.isfinite(pv.w_lor, "w_lor");
  hope(check_isfinite(pv.vel), "vel finite");
  hope(check_isfinite(pv.E), "E finite");
  hope(check_isfinite(pv.B), "B finite");
  
  return hope;
}

bool test_con2prim_mhd::check_isfinite(
          const cons_vars_mhd& cv,
          const prim_vars_mhd& pv) const
{
  return check_isfinite(pv) && check_isfinite(cv);
}


bool test_con2prim_mhd::chk_normal(real_t rho, real_t eps, real_t ye, 
                real_t z, real_t b, int vdim, int bdim) const
{
  failcount hope("Normal C2P Operation");

  prim_vars_mhd pv0;
  cons_vars_mhd cv0;
  setup_prim_cons(pv0, cv0, rho, eps, ye, z, b, vdim, bdim);

  prim_vars_mhd pv1;
  cons_vars_mhd cv1 = cv0;
  con2prim_mhd::report rep;
  hope(check_c2p(pv1, cv1, rep), "Calling C2P succeeds");

  if (hope(!rep.failed(), "C2P reports success")) {
    hope(compare_prims(pv0, pv1), 
         "C2P result close to original primitives");
    hope(check_same(cv0,cv1),
         "C2P not changing conservatives");
    hope(!rep.adjust_cons, 
         "C2P report adjust const flag not set");
  } 
  else {
    hope.postmortem(rep.debug_message());
  }
  hope(!rep.set_atmo, 
       "C2P report set atmo flag false for density above atmo cut");

  if (!hope) {
    hope.postmortem(str(format(
          "rho=%.15e, eps=%.15e, ye=%.15e, "
          "z(%d)=%.15e, b(%d)=%.15e")
        % pv0.rho % pv0.eps % pv0.ye % vdim % z % bdim % b));
  }

  return hope;
}
 


bool test_con2prim_mhd::chk_eps_adj(real_t rho, real_t deps, real_t ye, 
          real_t z, real_t b, int vdim, int bdim) const
{
  failcount hope("C2P works with invalid (too low) energy");  
  
  prim_vars_mhd pv0;
  cons_vars_mhd cv0;
  eos_thermal::range rgeps = eos.range_eps(rho, ye);
  setup_prim_cons(pv0, cv0, rho, rgeps.min(), ye, z, b, vdim, bdim);

  cons_vars_mhd cv1 = cv0;
  real_t dtau   = cv1.dens * deps / pv0.w_lor;
  cv1.tau      -= dtau;

  con2prim_mhd::report rep;
  prim_vars_mhd pv1;
  hope(check_c2p(pv1, cv1, rep), "Calling C2P succeeds");

  
  if (hope(!rep.failed(), "C2P adjustment succeeds")) {
    hope(compare_prims(pv1, pv0), 
         "C2P adjustment of primitives as expected");
    hope(compare_cons(cv1, cv0, pow(v_from_z(z),2)),
         "C2P adjustment of conserved vars as expected");
  }
  else {
    hope.postmortem(rep.debug_message());
  }
  hope(rep.adjust_cons, 
       "C2P report adjust cons flag set after adjustment or failure");
  hope(!rep.set_atmo,
       "C2P report set atmo flag not set after adjustment");

  if (!hope) {
    hope.postmortem(str(format(
          "tau=%.15e - %.15e, rho=%.15e, eps0=%.15e, "
          "ye=%.15e, z(%d)=%.15e, b(%d)=%.15e")
        % cv0.tau % dtau % pv0.rho % pv0.eps % pv0.ye
        % vdim % z % bdim % b));
  }
  
  return hope;
}


bool test_con2prim_mhd::chk_large_eps_adj(real_t rho, bool strict, 
          real_t deps, real_t ye, real_t z, real_t b, 
          int vdim, int bdim) const
{
  failcount hope("C2P works with invalid (too large) energy");  
  
  prim_vars_mhd pv0;
  cons_vars_mhd cv0;
  eos_thermal::range rgeps = eos.range_eps(rho, ye);
  setup_prim_cons(pv0, cv0, rho, rgeps.max(), ye, z, b, vdim, bdim);

  cons_vars_mhd cv1 = cv0;
  real_t dtau   = cv1.dens * deps / pv0.w_lor;
  cv1.tau      += dtau;

  con2prim_mhd::report rep;
  prim_vars_mhd pv1;
  hope(check_c2p(pv1, cv1, rep), "Calling C2P succeeds");
  
  if (strict) {
    hope(rep.failed(), 
         "C2P correctly reports fail in strict regime");
  } else {
    hope(!rep.failed(), 
         "C2P allows adjusting in non-strict regime");
    if (rep.failed()) hope.postmortem(rep.debug_message());
  }
  
  hope(rep.adjust_cons, 
       "C2P adjustment implies report flag adjust cons");
  hope(!rep.set_atmo, 
       "C2P adjustment implies report flag set atmo false");

  if (!hope) {
    hope.postmortem(str(format(
            "tau=%.15e - %.15e, rho=%.15e, eps0=%.15e, "
            "ye=%.15e, z(%d)=%.15e, b(%d)=%.15e")
          % cv0.tau % dtau % pv0.rho % pv0.eps % pv0.ye 
          % vdim % z % bdim % b));
  }


  return hope;
}


bool test_con2prim_mhd::chk_atmo(real_t rho_fac, real_t eps, real_t ye, 
                real_t z, real_t b, int vdim, int bdim) const
{
  failcount hope("C2P deals with atmosphere cut");
    
  prim_vars_mhd pv0, pv1;
  cons_vars_mhd cv0, cv1;
  con2prim_mhd::report rep;

  setup_prim_cons(pv0, cv0, rho_fac * atmo.rho_cut, 
                  eps, ye, z, b, vdim, bdim);
  cv1 = cv0;

  hope(check_c2p(pv1, cv1, rep), "Calling C2P succeeds");
  
  if (hope(!rep.failed(), "C2P not failing below atmo cut")) {
    hope(rep.set_atmo, 
         "C2P report set atmo flag correct");
    atmo.set(pv0, cv0, g);
    hope(check_same(pv1, pv0),
         "C2P atmosphere primitive vars as expected");
    hope(check_same(cv1, cv0),
         "C2P atmosphere conserved vars as expected");
  } else {
    hope.postmortem(rep.debug_message()); 
  }
  hope(rep.adjust_cons,
       "C2P report adjust const flag set for atmo case");
  
  if (!hope) {
    hope.postmortem(str(format(
        "rho=%.15e, eps0=%.15e, "
        "ye=%.15e, z(%d)=%.15e, b(%d)=%.15e")
      % pv0.rho % pv0.eps % pv0.ye % vdim % z % bdim % b));
  }
  
  return hope;
}

bool test_con2prim_mhd::chk_fail_rho(cons_vars_mhd cv) const
{
  failcount hope("C2P fails when density too high"); 

  con2prim_mhd::report rep;
  prim_vars_mhd pv;

  hope(check_c2p(pv, cv, rep), "Calling C2P succeeds");

  hope(rep.failed(), 
       "C2P reports fail when density invalid");
  hope(rep.status == con2prim_mhd::report::RANGE_RHO,
       "C2P reports correct reason for fail is density");

  return hope;
}

bool test_con2prim_mhd::chk_fail_eps(cons_vars_mhd cv) const
{
  failcount hope("C2P fails when density too high"); 

  con2prim_mhd::report rep;
  prim_vars_mhd pv;
  
  hope(check_c2p(pv, cv, rep), "Calling C2P succeeds");
  
  hope(rep.failed(), 
       "C2P reports fail when energy invalid");
  hope(rep.status == con2prim_mhd::report::RANGE_EPS,
       "C2P reports correct reason for fail is energy");

  return hope;
}





test_con2prim_mhd make_env(const env_idealgas& e) 
{
  
  auto eos = make_eos_idealgas(1.0, e.eps_max, e.rho_max);

  const real_t atmo_ye{ 0.5 };
  const real_t atmo_cut{ e.atmo_rho * 1.01 };
  const real_t atmo_p{ 
    eos.at_rho_eps_ye(e.atmo_rho, e.atmo_eps, atmo_ye).press() 
  };

  atmosphere atmo{e.atmo_rho, e.atmo_eps, atmo_ye, atmo_p, atmo_cut};

  const real_t rho_strict{ e.c2p_strict * e.atmo_rho };
  const bool   ye_lenient{ false };
  const int    max_iter{ 30 };
  const real_t max_b{ 10. };
  con2prim_mhd cv2pv(eos, rho_strict, ye_lenient, e.c2p_zmax, max_b, 
                     atmo, e.c2p_acc, max_iter);

  sm_metric3 g;
  g.minkowski();
  
  return test_con2prim_mhd(eos, g, atmo, cv2pv);
};

test_con2prim_mhd make_env(const env_hybrideos& e) 
{    
  auto eos = load_eos_thermal(PATH_EOS_HYB, units::geom_solar());

  const real_t atmo_cut{ e.atmo_rho * 1.01 };
  const real_t atmo_ye{ 0.25 };
  const real_t atmo_eps{ eos.range_eps(e.atmo_rho, atmo_ye).min() };
  BOOST_REQUIRE(eos.is_rho_valid(e.atmo_rho));
  BOOST_REQUIRE(eos.is_rho_valid(atmo_cut));
  BOOST_REQUIRE(eos.is_ye_valid(atmo_ye));

  const real_t atmo_p{ 
    eos.at_rho_eps_ye(e.atmo_rho, atmo_eps, atmo_ye).press() 
  };
  
  atmosphere atmo{e.atmo_rho, atmo_eps, atmo_ye, atmo_p, atmo_cut};

  const real_t rho_strict{ e.c2p_strict * e.atmo_rho };
  const bool   ye_lenient{ false };
  const int    max_iter{ 30 };
  const real_t max_b{ 10. };
  con2prim_mhd cv2pv(eos, rho_strict, ye_lenient, e.c2p_zmax, max_b, 
                     atmo, e.c2p_acc, max_iter);  
  
  sm_metric3 g;
  g.minkowski();
  
  return test_con2prim_mhd(eos, g, atmo, cv2pv);
};


BOOST_AUTO_TEST_CASE( c2p_mhd_basic )
{
  failcount hope{"C2P basic tests"};
  
  env_idealgas par{1e6, 100., 1e-11, 1e-6, 10.0, 10.0, 4e-8};

  const auto tst = make_env(par);
  
  const real_t ye_fixed{0.25};
  const real_t rho{1e-5};
  const real_t rho_nostrict{ 
     par.atmo_rho * ( 0.1 + 0.9*par.c2p_strict)}; 

  hope(tst.chk_normal(rho, 0.01, ye_fixed, 0.01, 0, 2, 1),
       "C2P works for eps and v small, B=0");

  hope(tst.chk_normal(rho, 10.0, ye_fixed, 0.01, 0, 2, 1),
       "C2P works for eps large, v small, B=0");

  hope(tst.chk_normal(rho, 10.0, ye_fixed, 9.9, 0, 2, 1),
       "C2P works for eps large, v large, B=0");

  hope(tst.chk_normal(rho, 0.01, ye_fixed, 9.9, 0, 2, 1),
       "C2P works for eps small, v large, B=0");
             
  hope(tst.chk_normal(rho, 0.01, ye_fixed, 0.01, 5, 2, 1),
       "C2P works for eps small, v small, B large");
                  
  hope(tst.chk_normal(rho, 10, ye_fixed, 0.01, 5, 2, 1),
       "C2P works for eps large, v small, B large");

  hope(tst.chk_normal(rho, 10, ye_fixed, 9.9, 5, 2, 1),
       "C2P works for eps large, v large, B large");
                  
  hope(tst.chk_normal(rho, 0.316228, ye_fixed, 
                      1., 3.5, 2, 1),
       "C2P works for moderate eps, v, B");
       
  hope(tst.chk_eps_adj(rho, 0.1, ye_fixed, 1., 3.5, 2, 1),
       "C2P adjusts epsilon that is moderately too small");
          
         
  hope(tst.chk_normal(rho_nostrict, 99, ye_fixed, 1., 3.5, 2, 1),
       "C2P works in non-strict regime for huge eps moderate v");

  hope(tst.chk_large_eps_adj(rho_nostrict, false, 2.0, 
                         ye_fixed, 1., 3.5, 2, 1),
       "C2P works for eps above EOS range in non-strict regime");
              
  hope(tst.chk_large_eps_adj(rho, true, 2.0, ye_fixed, 
                             1., 3.5, 2, 1),
       "C2P works for eps above EOS range in strict regime");

  hope(tst.chk_atmo(0.9, 5.0, ye_fixed, 0.01, 3.5, 2, 1),
       "C2P works for rho and dens below atmo cut");            
  
  hope(tst.chk_atmo(0.9, 5.0, ye_fixed, 9., 3.5, 2, 1),
       "C2P works for rho below atmo cut, but dens above");
}


BOOST_AUTO_TEST_CASE( c2p_mhd_fail )
{
  failcount hope{"C2P really fails when input violates error policy"}; 
  
  env_idealgas par{1e-3, 10., 1e-12, 1e-6, 10.0, 10.0, 4e-8};
  
  const auto tst = make_env(par);

  cons_vars_mhd cv1{2e-3, 1e-4, 1e-3, {0.,0.,0.}, {0.,0.,0.}};    
  hope(tst.chk_fail_rho(cv1),
       "C2P fails if dens much too high (at zero velocity)");
          
  cons_vars_mhd cv2{1e-6, 2e-5, 5e-7, {0.,0.,0.}, {0.,0.,0.}};   
  hope(tst.chk_fail_eps(cv2),
       "C2P fails if tau much too high (at zero velocity)");
  
}


BOOST_AUTO_TEST_CASE( test_con2prim_phys_igas )
{
  failcount hope{"C2P with ideal gas EOS works "
                 "in large parameter space"};
    
  env_idealgas par{1e6, 51., 1e-11, 1e-6, 1.0, 2e3, 1e-8};
  
  const auto tst = make_env(par);  
                 
  const real_t max_z{ 1e3 };
  const real_t min_z{ 1e-2 };
  const real_t max_b{ 5.0 };
  const auto iz = log_spacing(min_z, max_z, 100);
  const auto ieps = linear_spacing(1E-4, 50., 10);
  const auto ib = linear_spacing(0.0, max_b, 10);

  const real_t ye_fixed{0.25};
  const real_t rho_fixed{1e-5};

  for (int vdim=0; vdim<3; vdim++) {
    for (int bdim=0; bdim<3; bdim++) {
      for (const real_t z : iz) {
        for (const real_t eps : ieps) {
          for (const real_t b : ib) {
            hope(tst.chk_normal(rho_fixed, eps, ye_fixed, z, b, 
                                vdim, bdim), 
                 "C2P works at valid point in parameter space");
          }
        }
      }
    }
  }

}

BOOST_AUTO_TEST_CASE( test_con2prim_phys_hybr )
{
  failcount hope("C2P with hybrid EOS works "
                 "in large parameter space");

  env_hybrideos par{1e-12, 1.0, 2e3, 1e-8};
  
  const auto tst = make_env(par);             
  
  real_t max_z{ 1e3};
  real_t min_z{ 1e-2};
  real_t max_b{ 5.0 };

  const real_t ye_fixed{0.25};

  auto iz   = log_spacing(min_z, max_z, 100);
  auto ith  = linear_spacing(2e-6, 0.99, 10); //2e-6 = 1e-4 / epsmax
  auto irho = log_spacing(par.atmo_rho*5, 
                          tst.eos.range_rho().max()/1.1, 10);
  auto ib   = linear_spacing(0.0, max_b, 10);
  
  int vdim=0;
  const real_t ye = ye_fixed;
  for (const real_t rho : irho) {
    eos_thermal::range rgeps = tst.eos.range_eps(rho, ye);
    for (const real_t z : iz) {
      for (const real_t w : ith) {          
        const real_t eps = rgeps.min()*(1.0-w) + rgeps.max()*w;
        for (int bdim=0; bdim<2; bdim++) {
          for (const real_t b : ib) {
            hope(tst.chk_normal(rho, eps, ye, z, b, vdim, bdim),
                 "C2P works at valid point in parameter space");
          }
        }
      }
    }
  }
  //~ }
}



