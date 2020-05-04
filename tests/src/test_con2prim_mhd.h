#include "test_utils.h"
#include "con2prim_imhd.h"
#include "eos_thermal.h"
#include <boost/format.hpp>

using EOS_Toolkit::eos_thermal;
using EOS_Toolkit::real_t;
using EOS_Toolkit::sm_metric3;
using EOS_Toolkit::sm_tensor1;
using EOS_Toolkit::atmosphere;
using EOS_Toolkit::con2prim_mhd;
using EOS_Toolkit::prim_vars_mhd;
using EOS_Toolkit::cons_vars_mhd;



struct env_idealgas {
  real_t rho_max; 
  real_t eps_max; 
  real_t atmo_rho; 
  real_t atmo_eps;
  real_t c2p_strict; 
  real_t c2p_zmax;
  real_t c2p_acc;
};


struct env_hybrideos {
  real_t atmo_rho; 
  real_t c2p_strict; 
  real_t c2p_zmax;
  real_t c2p_acc;
};

class test_con2prim_mhd {
  public:

  eos_thermal eos;
  sm_metric3 g;
  atmosphere atmo;
  con2prim_mhd cv2pv;
  
  
  test_con2prim_mhd(const eos_thermal& eos_,
                    const sm_metric3& g_,
                    const atmosphere &atmo_,
                    const con2prim_mhd& cv2pv_)
  : eos(eos_), g(g_), atmo(atmo_), cv2pv(cv2pv_) {}
                    
  
  void setup_prim_cons(prim_vars_mhd& pv, 
                       cons_vars_mhd& cv, 
                       real_t rho, real_t eps, real_t ye, 
                       real_t z, real_t b, int vdim, int bdim) const;

  bool check_c2p(prim_vars_mhd& pv, 
                 cons_vars_mhd& cv, 
                 con2prim_mhd::report& rep) const;

  bool compare_prims(const prim_vars_mhd&pv0, 
                     const prim_vars_mhd& pv1) const;

  bool compare_cons(const cons_vars_mhd& cv0, 
                    const cons_vars_mhd& cv1, real_t vsqr) const;

  template<class T, int N, bool UP> 
  bool check_isclose(const sm_tensor1<T,N,UP>& v1,
                     const sm_tensor1<T,N,UP>& v2, 
                     real_t reltol, real_t abstol) const;

  bool check_same(const prim_vars_mhd&pv0, 
                  const prim_vars_mhd& pv1) const;

  bool check_same(const cons_vars_mhd& cv0, 
                  const cons_vars_mhd& cv1) const;

  bool check_vel_w(const prim_vars_mhd&pv) const;
                   
  template<class T, int N, bool UP> 
  bool check_isnan(const sm_tensor1<T,N,UP>& v) const;
                      
  bool check_isnan(const cons_vars_mhd& cv) const;
  bool check_isnan(const prim_vars_mhd& pv) const;
  bool check_isnan(const cons_vars_mhd& cv,
                   const prim_vars_mhd& pv) const;

  template<class T, int N, bool UP> 
  bool check_isfinite(const sm_tensor1<T,N,UP>& v) const;

  bool check_isfinite(const cons_vars_mhd& cv) const;
  bool check_isfinite(const prim_vars_mhd& pv) const;
  bool check_isfinite(const cons_vars_mhd& cv,
                      const prim_vars_mhd& pv) const;


  bool chk_normal(real_t rho, real_t eps, real_t ye, 
                  real_t z, real_t b, int vdim, int bdim) const;
                  
  bool chk_fail_rho(cons_vars_mhd cv) const;
  bool chk_fail_eps(cons_vars_mhd cv) const;
                  

  bool chk_eps_adj(real_t rho, real_t deps, real_t ye, 
                   real_t z, real_t b, int vdim, int bdim) const;

  bool chk_large_eps_adj(real_t rho, bool strict, real_t deps, 
                         real_t ye, real_t z, real_t b, 
                         int vdim, int bdim) const;
  
  bool chk_atmo(real_t rho_fac, real_t eps, real_t ye, 
                real_t z, real_t b, int vdim, int bdim) const;
};


template<class T, int N, bool UP> 
bool test_con2prim_mhd::check_isnan(const sm_tensor1<T,N,UP>& v) const
{
  failcount hope("vector components all NAN");
  
  for (int d=0; d<N; d++) {
    hope.isnan(v(d), (boost::format("v(%d)") % d).str());  
  }
  return hope;
}

template<class T, int N, bool UP> 
bool 
test_con2prim_mhd::check_isfinite(const sm_tensor1<T,N,UP>& v) const
{
  failcount hope("vector components all finite");
  
  for (int d=0; d<N; d++) {
    hope.isfinite(v(d), (boost::format("v(%d)") % d).str());  
  }
  return hope;
}

template<class T, int N, bool UP> 
bool test_con2prim_mhd::check_isclose(
             const sm_tensor1<T,N,UP>& v1,
             const sm_tensor1<T,N,UP>& v2, 
             real_t reltol, real_t abstol) const
{
  failcount hope("Vectors close within tolerance");
  
  const real_t r{ sqrt(std::max(g.norm2(v1), g.norm2(v2))) };
  const real_t dv{ sqrt(g.norm2(v1-v2)) };
  if (!hope.isleq(dv, r*reltol+abstol, "norm of difference vector")) {
    for (int d=0; d<N; d++) {
      hope.postmortem((boost::format("v1(%1%)=%2$.15e, v2(%1%)=%3$.15e") 
                       % d % v1(d) % v2(d)).str());  
    }    
  }
  
  return hope;
}

real_t w_from_z(real_t z);
real_t v_from_z(real_t z);
