#include "c2p_report.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <stdexcept>
#include <cmath>

using namespace std;

namespace EOS_Toolkit {


void c2p_mhd_report::raise() const
{
  throw runtime_error(debug_message()); 
}

void c2p_mhd_report::set_atmo_set()
{
  status      = SUCCESS;
  set_atmo    = true;
  adjust_cons = true;
}

void c2p_mhd_report::set_invalid_detg(real_t detg_)
{
  status      = INVALID_DETG;
  adjust_cons = true;
  detg        = detg_;
}

void c2p_mhd_report::set_nans_in_cons(real_t dens_, real_t qtot_, 
  real_t rsqr_, real_t rbsqr_, real_t bsqr_, real_t ye_)
{
  status      = NANS_IN_CONS;
  set_atmo    = false;
  adjust_cons = true;
  dens        = dens_;
  qtot        = qtot_;
  rsqr        = rsqr_;
  rbsqr       = rbsqr_;
  bsqr        = bsqr_;
  ye          = ye_;
}

void c2p_mhd_report::set_range_rho(real_t dens_, real_t rho_)
{
  status      = RANGE_RHO;
  set_atmo    = false;
  adjust_cons = true;
  dens        = dens_;
  rho         = rho_;
}
  
void c2p_mhd_report::set_range_eps(real_t eps_)
{
  status      = RANGE_EPS;
  set_atmo    = false;
  adjust_cons = true;
  eps         = eps_;
}

void c2p_mhd_report::set_speed_limit(real_t vel_)
{
  status      = SPEED_LIMIT;
  set_atmo    = false;
  adjust_cons = true;
  vel         = vel_;
}

void c2p_mhd_report::set_b_limit(real_t bsqr_)
{
  status      = B_LIMIT;
  set_atmo    = false;
  adjust_cons = true;
  bsqr        = bsqr_;
}

void c2p_mhd_report::set_range_ye(real_t ye_)
{
  status      = RANGE_YE;
  set_atmo    = false;
  adjust_cons = true;
  ye          = ye_;
}

void c2p_mhd_report::set_root_conv()
{
  status      = ROOT_FAIL_CONV;
  set_atmo    = false;
  adjust_cons = true;
}

void c2p_mhd_report::set_root_bracket()
{
  status      = ROOT_FAIL_BRACKET;
  set_atmo    = false;
  adjust_cons = true;
}

void c2p_mhd_report::set_prep_root_conv()
{
  status          = PREP_ROOT_FAIL_CONV;
  set_atmo        = false;
  adjust_cons     = true;
}

void c2p_mhd_report::set_prep_root_bracket()
{
  status          = PREP_ROOT_FAIL_BRACKET;
  set_atmo        = false;
  adjust_cons     = true;
}


string c2p_mhd_report::debug_message() const 
{
  ostringstream os;
  os << scientific << setprecision(15);
  switch (status) {
    case SUCCESS:
      os << "Con2Prim succeeded.";
      if (set_atmo) {
        os << "Artificial atmosphere has been enforced.";
      }
      if (adjust_cons) {
        os << "Conserved variables have been changed.";
      }
      break;
    case INVALID_DETG:
      os << "Invalid metric determinant (" << detg << ")";
      break;
    case NANS_IN_CONS:
      os << "NAN in conserved variables"
         << ", dens=" << dens 
         << ", qtot=" << qtot 
         << ", rsqr=" << rsqr 
         << ", rbsq=" << rbsqr 
         << ", bsqr=" << bsqr 
         << ", ye=" << ye;
      break;
    case RANGE_RHO: 
      os << "Density out of range, dens=" << dens
         << ", rho=" << rho;
      break;
    case RANGE_EPS:
      os << "Specific energy out of range, eps=" << eps;
      break;
    case SPEED_LIMIT:
      os << "Speed limit exceeded, v=" << vel;
      break;
    case B_LIMIT:
      os << "Limit for magnetic field exceeded, b=" << std::sqrt(bsqr);
      break;
    case RANGE_YE:
      os << "Electron fraction out of range, Y_e=" << ye;
      break;
    case ROOT_FAIL_CONV:
      os << "Root finding failed (not converged after " 
         << iters << "steps)";
      break;
    case ROOT_FAIL_BRACKET:
      os << "Root finding failed (faulty bracketing)";
      break;
    case PREP_ROOT_FAIL_CONV:
      os << "Preparatory root finding failed (not converged)";
      break;
    case PREP_ROOT_FAIL_BRACKET:
      os << "Preparatory root finding failed (faulty bracketing)";
      break;
    default:
      assert(false);
      os << "Invalid error type. Should never happen. Code is messed up.";
      break;
  }
  return os.str();
}

}

