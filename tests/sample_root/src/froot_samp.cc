#include "sroot_config.h"
#include <cassert>
#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "unitconv.h"
#include "con2prim_imhd_internals.h"
#include "eos_thermal.h"
#include "eos_idealgas.h"
#include "eos_thermal_file.h"
#include "eos_hybrid.h"

using namespace std;
using namespace EOS_Toolkit;


struct prims {
  string name;
  real_t rho;
  real_t eps_th;
  real_t ye;
  array<real_t,3> vel;
  real_t b;
};


void setup_cons(cons_vars_mhd& cv, prim_vars_mhd &pv,
                     const sm_metric3& g, 
                     const eos_thermal& eos, 
                     real_t rho, real_t eps_th, real_t ye, 
                     array<real_t,3> vel, real_t b)
{
  real_t eps = eps_th + eos.range_eps(rho, ye).min();
  
  real_t press = eos.at_rho_eps_ye(rho, eps, ye).press();
  
  sm_vec3u v;
  for (int i =0; i<3; i++) v(i) = vel[i];
  real_t vsqr = g.contract(v,v);
  real_t w    = 1.0 / sqrt(1.0 - vsqr);
  sm_vec3u B   = ZERO;
  B(0)      = b*sqrt(rho)*w;
  sm_vec3l El = g.cross_product(B, v);
  sm_vec3u E   = g.raise(El);
  pv            = prim_vars_mhd(rho, eps, ye, press, v, w, E, B);
  cv.from_prim(pv, g);
}


void sample_froot(const eos_thermal eos,   
                   const prims& p, string path, int nsamp=10000) {
  sm_metric3 g;
  g.minkowski();
  cons_vars_mhd cv;
  prim_vars_mhd pv;
  setup_cons(cv, pv, g, eos, p.rho, p.eps_th, p.ye, p.vel, p.b);  
  const real_t h = 1 + pv.eps + pv.press / pv.rho;
  const real_t mu_0 = 1.0 / (pv.w_lor * h);
  
  const real_t d     = cv.dens / g.vol_elem;
  const sm_vec3u bu   = cv.bcons / (g.vol_elem * sqrt(d));
  const sm_vec3l rl   = cv.scon / cv.dens;

  const sm_vec3u ru   = g.raise(rl);
  const real_t rsqr  = ru * rl;
  const real_t rb    = rl * bu;
  const real_t rbsqr = rb * rb;
  const real_t bsqr  = g.contract(bu, bu);
  const real_t q     = cv.tau / cv.dens;
  const real_t ye0   = cv.tracer_ye / cv.dens;

  
  detail::froot::cache sol;
  detail::froot f(eos, ye0, d, q, rsqr, rbsqr, bsqr, sol); 

  con2prim_mhd::report errs;
  auto bracket = f.initial_bracket(errs);  
  
  ofstream os(path + "/" + p.name.c_str() + ".dat"); 
  os << setprecision(15);
  
  for (int i=0; i<=nsamp; i++) {
    const real_t mu = bracket.min() 
                      + bracket.length() * real_t(i) / nsamp; 
    const real_t dmu = mu - mu_0;
    os << setw(30) <<left 
       << mu << " " 
       << dmu << " " 
       << f(mu) << endl;
  }
  
}

int main(int argc, char *argv[])
{
  assert(argc==2);

  
  string eos_cold_file = PATH_EOS_HYB;
  
  auto eos = load_eos_thermal(PATH_EOS_HYB, units::geom_solar());
  const real_t rho_max = eos.range_rho().max();
   
  
  array<prims,8> example_prims{{
    {"cold_vzero_bzero", rho_max/5, 0, 0.5, {0.,0.,0.}, 0.},
    {"cold_vlarge_bzero", rho_max/5, 0, 0.5, {0., 0.99, 0.}, 0.},
    {"cold_vzero_blarge", rho_max/5, 0, 0.5, {0.,0.,0.}, 1.},
    {"cold_vlarge_blarge", rho_max/5, 0, 0.5, {0., 0.99, 0.}, 1.},
    {"hot_vzero_bzero", rho_max/5, 10., 0.5, {0.,0.,0.}, 0.},
    {"hot_vzero_blarge", rho_max/5, 10., 0.5, {0.,0.,0.}, 1.},
    {"hot_vlarge_bzero", rho_max/5, 10., 0.5, {0., 0.99, 0.}, 0.},
    {"hot_vlarge_blarge", rho_max/5, 10., 0.5, {0., 0.99, 0.}, 1.}
  }};



  for (const prims& p : example_prims) {
    sample_froot(eos, p, argv[1]);
  }

  return 0;
}


