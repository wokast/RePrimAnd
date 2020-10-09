#include "bench_config.h"
#include "bench_utils.h"

#include <cassert>
#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "con2prim_imhd.h"
#include "eos_thermal_file.h"
#include "eos_hybrid.h"
#include "eos_idealgas.h"

using namespace std;
using namespace EOS_Toolkit;


real_t w_from_z(real_t z) {
  return sqrt(1.0+z*z);
}

real_t v_from_z(real_t z) {
  return z/w_from_z(z);
}

void setup_prim_cons(EOS_Toolkit::prim_vars_mhd& pv, 
                     EOS_Toolkit::cons_vars_mhd& cv, 
                     const EOS_Toolkit::sm_metric3& g, 
                     const EOS_Toolkit::eos_thermal& eos, 
                     real_t rho, real_t eps, real_t ye, real_t v, 
                     real_t b, int vdim, int bdim)
{
  real_t press0 = 0.0;
  if (rho > 0) {
    press0 = eos.at_rho_eps_ye(rho, eps, ye).press();
  }
  sm_vec3u vel0 = ZERO;
  vel0(vdim)    = v;
  real_t wl0    = 1.0 / sqrt(1.0 - v*v);
  sm_vec3u B0   = ZERO;
  B0(bdim)      = b*sqrt(rho*wl0);
  sm_vec3l El   = g.cross_product(B0, vel0);
  sm_vec3u E0   = g.raise(El);

  pv            = prim_vars_mhd(rho, eps, ye, press0, vel0, wl0, E0, B0);
  cv.from_prim(pv, g);
}




void map_z_eps(eos_thermal eos, con2prim_mhd& cv2pv, real_t b, real_t rho, 
               real_t ye, sm_metric3& g, string fn) 
{
  auto iz     = log_spacing(1e-2, 1e3, 200);
  auto ieps   = log_spacing(1e-4, 1e1, 200);
  const real_t eps0     = eos.range_eps(rho, ye).min();
  ofstream os(fn.c_str()); 
  
  for (const real_t z : iz) {
    const real_t v = v_from_z(z);
    for (const real_t epsth : ieps) {
      const real_t eps = eps0 + epsth;

      prim_vars_mhd pv;
      cons_vars_mhd cv;
      setup_prim_cons(pv, cv, g, eos, rho, eps, ye, v, b, 0, 1);
      con2prim_mhd::report rep;
      cv2pv(pv, cv, g, rep);
      if (rep.failed()) {
        cout << rep.debug_message() <<endl;
        assert(false);
      }
      assert(!rep.failed());
      os << setw(20) << z << setw(20) << epsth 
         << setw(20) << rep.iters 
         << setw(20) << pv.press 
         << setw(20) << b*sqrt(cv.dens)
         << endl;
    }
    os << endl;
  }
}

void map_z_b(eos_thermal eos, con2prim_mhd& cv2pv, real_t epsth, 
             real_t rho, real_t ye, sm_metric3& g, string fn) 
{
  auto iz     = log_spacing(1e-2, 1e3, 200);
  auto ib     = log_spacing(1e-4, 1e4, 200);
  
  const real_t eps0  = eos.range_eps(rho, ye).min();
  const real_t eps   = eps0 + epsth;
  ofstream os(fn.c_str());  

  for (const real_t z : iz) {
    const real_t v = v_from_z(z);
    for (const real_t b : ib) {
      prim_vars_mhd pv;
      cons_vars_mhd cv;
      setup_prim_cons(pv, cv, g, eos, rho, eps, ye, v, b, 0, 1);
      con2prim_mhd::report rep;
      cv2pv(pv, cv, g, rep);
      if (rep.failed()) {
        cout << rep.debug_message() <<endl;
        assert(false);
      }
      assert(!rep.failed());
      os << setw(20) << z << setw(20) << b <<setw(20) << rep.iters 
         << setw(20) << pv.press 
         << setw(20) << b*sqrt(cv.dens)
         << endl;
    }
    os << endl;
  }
}

atmosphere get_atmo(eos_thermal eos, const real_t eps_th=0.,
                    const real_t rho_atmo= 1e-11,
                    const real_t ye_atmo= 0.25) 
{
  const real_t rho_atmo_cut = rho_atmo * 1.01;
  assert(eos.is_rho_valid(rho_atmo));
  assert(eos.is_rho_valid(rho_atmo_cut));
  assert(eos.is_ye_valid(ye_atmo));
  
  const real_t eps0      = eos.range_eps(rho_atmo, ye_atmo).min();
  const real_t eps_atmo  = eps_th + eps0;
  
  const real_t p_atmo       = eos.at_rho_eps_ye(rho_atmo, 
                                         eps_atmo, ye_atmo). press();
  
  return atmosphere(rho_atmo, eps_atmo, ye_atmo, p_atmo, rho_atmo_cut);
}


eos_thermal get_eos_ig() {
  real_t max_eps = 100.;
  real_t max_rho = 1e6;
  return make_eos_idealgas(1.0, max_eps, max_rho);
}


int main(int argc, char *argv[])
{
  assert(argc==2);
  const string path{string(argv[1])+"/"};
  
  eos_thermal eos_ig = get_eos_ig();
  atmosphere atmo_ig = get_atmo(eos_ig, 1e-6);
  const real_t acc          = 1e-8;  
  con2prim_mhd cv2pv_ig(eos_ig, atmo_ig.rho, false, 2e3, 5e4, 
                        atmo_ig, acc, 100);


  auto eos_hyb = load_eos_thermal(PATH_EOS_HYB, units::geom_solar());

  atmosphere atmo_hyb = get_atmo(eos_hyb, 0.0);
  con2prim_mhd cv2pv_hyb(eos_hyb, atmo_hyb.rho, false, 2e3, 5e4, 
                         atmo_hyb, acc, 100);


  sm_metric3 g;
  g.minkowski();

  const real_t ye_fixed  = 0.25;
  const real_t rho_fixed = 1e-5;  
  const real_t blarge    = 10.;
  const real_t epsth_hot = 10.;
  const real_t epsth_cold = 1e-4;
  
  map_z_eps(eos_ig, cv2pv_ig, 0, rho_fixed, ye_fixed, g, 
            path+"perf_eosig_z_eps_Bzero.dat");
  map_z_eps(eos_ig, cv2pv_ig, blarge, rho_fixed, ye_fixed, g, 
            path+"perf_eosig_z_eps_Blarge.dat");
  map_z_b(eos_ig, cv2pv_ig, epsth_cold, rho_fixed, ye_fixed, g, 
            path+"perf_eosig_z_b_cold.dat");
  map_z_b(eos_ig, cv2pv_ig, epsth_hot, rho_fixed, ye_fixed, g, 
            path+"perf_eosig_z_b_hot.dat");


  map_z_eps(eos_hyb, cv2pv_hyb, 0, rho_fixed, ye_fixed, g, 
            path+"perf_eoshyb_z_eps_Bzero.dat");
  map_z_eps(eos_hyb, cv2pv_hyb, blarge, rho_fixed, ye_fixed, g, 
            path+"perf_eoshyb_z_eps_Blarge.dat");
  map_z_b(eos_hyb, cv2pv_hyb, epsth_cold, rho_fixed, ye_fixed, g, 
            path+"perf_eoshyb_z_b_cold.dat");
  map_z_b(eos_hyb, cv2pv_hyb, epsth_hot, rho_fixed, ye_fixed, g, 
            path+"perf_eoshyb_z_b_hot.dat"); 
            
  return 0;
}


