#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include <stdexcept>
#include <limits>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "reprimand/eos_thermal.h"
#include "reprimand_global_eos.h"

using namespace EOS_Toolkit;
using namespace RePrimAnd_Global_EOS;


namespace { 
  
const CCTK_INT eos_key_reprimand_thermal = 13666;
const CCTK_INT eos_key_reprimand_cold    = 13766;


void warn_invalid_handle(CCTK_INT eoskey)
{
  CCTK_WARN(0, "Invalid EOS handle in RePrimAnd_EOS_Omni_API");
}


void warn_invalid_rho_eps_ye(const eos_thermal& eos, 
  CCTK_REAL rho, CCTK_REAL eps, CCTK_REAL ye)
{  
  DECLARE_CCTK_PARAMETERS;
  if (0 == verbosity) return;

  std::stringstream msg;
  msg.precision(16);
  msg << "RePrimAnd_EOS_Omni_API: invalid input" << std::endl    
      << "rho = " << rho << ", "
      << "eps = " << eps << ", "
      << "Y_e = " << ye << std::endl;
  auto rgrho = eos.range_rho();
  auto rgye = eos.range_ye();

  if (rgrho.contains(rho)) {
    if (rgye.contains(ye)) {
      auto rgeps = eos.range_eps(rho, ye);
      if (!rgeps.contains(eps)) {
        msg << "eps out of valid range [" << rgeps.min() << ", " 
            << rgeps.max() << "]" << std::endl;
        if (sloppy_eps) {
          msg << "Adjusted eps to valid range" << std::endl;
        }
      }
    }
  }
  else {
    msg << "rho out of valid range [" << rgrho.min() << ", " 
    << rgrho.max() << "]" << std::endl;
  }
  if (!rgye.contains(ye)) {
    msg << "Y_e out of valid range [" << rgye.min() << ", " 
        << rgye.max() << "]" << std::endl;
  }
  CCTK_WARN(1, msg.str().c_str());
}

eos_thermal::state state_rho_eps_ye(const eos_thermal& eos, 
  CCTK_REAL rho, CCTK_REAL eps, CCTK_REAL ye, CCTK_INT& errcode)
{
  DECLARE_CCTK_PARAMETERS;

  auto s = eos.at_rho_eps_ye(rho, eps, ye);
  if (s) {
    errcode = 0;
    return s;
  }

  warn_invalid_rho_eps_ye(eos, rho, eps, ye);

  if (sloppy_eps && eos.is_rho_ye_valid(rho, ye)) 
  {
    real_t valid_eps = eos.range_eps(rho, ye).limit_to(eps);
    
    errcode = 0; 
    return eos.at_rho_eps_ye(rho, valid_eps, ye);
  }
  
  errcode = -1;
  return s;
}


void warn_invalid_rho_temp_ye(const eos_thermal& eos, 
  CCTK_REAL rho, CCTK_REAL temp, CCTK_REAL ye)
{  
  DECLARE_CCTK_PARAMETERS;
  if (0 == verbosity) return;

  std::stringstream msg;
  msg.precision(16);
  msg << "RePrimAnd_EOS_Omni_API: invalid input" << std::endl    
      << "rho = " << rho << ", "
      << "T   = " << temp << ", "
      << "Y_e = " << ye << std::endl;
  auto rgrho = eos.range_rho();
  auto rgye = eos.range_ye();

  if (rgrho.contains(rho)) {
    if (rgye.contains(ye)) {
      auto rgt = eos.range_temp(rho, ye);
      if (!rgt.contains(temp)) {
        msg << "T out of valid range [" << rgt.min() << ", " 
            << rgt.max() << "]" << std::endl;
        if (sloppy_temp) {
          msg << "Adjusted T to valid range" << std::endl;
        }
      }
    }
  }
  else {
    msg << "rho out of valid range [" << rgrho.min() << ", " 
    << rgrho.max() << "]" << std::endl;
  }
  if (!rgye.contains(ye)) {
    msg << "Y_e out of valid range [" << rgye.min() << ", " 
        << rgye.max() << "]" << std::endl;
  }
  CCTK_WARN(1, msg.str().c_str());
}

eos_thermal::state state_rho_temp_ye(const eos_thermal& eos, 
  CCTK_REAL rho, CCTK_REAL temp, CCTK_REAL ye, CCTK_INT& errcode)
{
  DECLARE_CCTK_PARAMETERS;

  auto s = eos.at_rho_temp_ye(rho, temp, ye);
  if (s) {
    errcode = 0;
    return s;
  }

  warn_invalid_rho_temp_ye(eos, rho, temp, ye);

  if (sloppy_temp && eos.is_rho_ye_valid(rho, ye)) 
  {
    real_t valid_temp = eos.range_temp(rho, ye).limit_to(temp);
    
    errcode = 0; 
    return eos.at_rho_temp_ye(rho, valid_temp, ye);
  }
  
  errcode = -1;
  return s;
}


void warn_invalid_rho(const eos_barotr& eos, CCTK_REAL rho)
{  
  DECLARE_CCTK_PARAMETERS;
  if (0 == verbosity) return;

  std::stringstream msg;
  msg.precision(16);
  msg << "RePrimAnd_EOS_Omni_API: invalid input" << std::endl    
      << "rho = " << rho << std::endl;
  auto rgrho = eos.range_rho();

  msg << "rho out of valid range [" << rgrho.min() << ", " 
      << rgrho.max() << "]" << std::endl;
  
  CCTK_WARN(1, msg.str().c_str());
}

eos_barotr::state state_rho(const eos_barotr& eos, 
                             CCTK_REAL rho, CCTK_INT& errcode)
{
  auto s = eos.at_rho(rho);
  
  if (s) {
    errcode = 0;
  } 
  else {
    warn_invalid_rho(eos, rho);
    errcode = -1;
  }
  
  return s;
}




void press_eps_from_rho(const CCTK_INT npoints, 
  const CCTK_REAL* rho, CCTK_REAL* eps, CCTK_REAL* press,
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{

  const eos_barotr& eos = global_eos_cold::get_eos();
  *anyerr = 0;
  for (int i=0; i<npoints; ++i) 
  {
    auto s = state_rho(eos, rho[i], keyerr[i]);
    if (s) 
    {
      eps[i]    = s.eps();
      press[i]  = s.press();
    }
    else 
    {
      eps[i]    = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      press[i]  = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      *anyerr   = 1;
    }
  }
}

void press_eps_from_rho_temp_ye(const CCTK_INT npoints, 
  const CCTK_REAL* rho, CCTK_REAL* eps, const CCTK_REAL* temp,
  const CCTK_REAL* ye, CCTK_REAL* press,
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  const eos_thermal& eos = global_eos_thermal::get_eos();
  *anyerr = 0;
  for (int i=0; i<npoints; ++i) 
  {
    auto s = state_rho_temp_ye(eos, rho[i], temp[i], ye[i], keyerr[i]);
    if (s) 
    {
      eps[i]    = s.eps();
      press[i]  = s.press();
    }
    else 
    {
      eps[i]    = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      press[i]  = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      *anyerr   = 1;
    }
  }    
}



void press_from_rho_eps_ye(const CCTK_INT npoints, 
  const CCTK_REAL* rho, const CCTK_REAL* eps, 
  const CCTK_REAL* ye, CCTK_REAL* press,
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{

  const eos_thermal& eos = global_eos_thermal::get_eos();
  *anyerr = 0;
  for (int i=0; i<npoints; ++i) 
  {
    auto s = state_rho_eps_ye(eos, rho[i], eps[i], ye[i], keyerr[i]);
    
    if (s) 
    {
      press[i]  = s.press();
    }
    else 
    {
      press[i]  = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      *anyerr   = 1;
    }
  }
}


void press_eps_cs2_from_rho(const CCTK_INT npoints, 
  const CCTK_REAL* rho, CCTK_REAL* eps, 
  CCTK_REAL* press, CCTK_REAL* cs2,
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  const eos_barotr& eos = global_eos_cold::get_eos();
  *anyerr = 0;
  for (int i=0; i<npoints; ++i) 
  {
    auto s = state_rho(eos, rho[i], keyerr[i]);
    if (s) 
    {
      eps[i]    = s.eps();
      press[i]  = s.press();
      CCTK_REAL cs = s.csnd();
      cs2[i]    = cs*cs;
    }
    else 
    {
      eps[i]    = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      press[i]  = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      cs2[i]    = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      *anyerr   = 1;
    }
  }
}



void press_eps_cs2_from_rho_temp_ye(const CCTK_INT npoints, 
  const CCTK_REAL* rho, CCTK_REAL* eps, const CCTK_REAL* temp, 
  const CCTK_REAL* ye, CCTK_REAL* press, CCTK_REAL* cs2,
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{  
  const eos_thermal& eos = global_eos_thermal::get_eos();
  *anyerr = 0;
  for (int i=0; i<npoints; ++i) 
  {
    auto s = state_rho_temp_ye(eos, rho[i], temp[i], ye[i], keyerr[i]);
    if (s) 
    {
      eps[i]    = s.eps();
      press[i]  = s.press();
      CCTK_REAL cs = s.csnd();
      cs2[i]    = cs*cs;
    }
    else 
    {
      eps[i]    = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      press[i]  = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      cs2[i]    = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      *anyerr   = 1;
    }
  }
}


void press_cs2_from_rho_eps_ye(const CCTK_INT npoints, 
  const CCTK_REAL* rho, const CCTK_REAL* eps, 
  const CCTK_REAL* ye, CCTK_REAL* press, CCTK_REAL* cs2,
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  const eos_thermal& eos = global_eos_thermal::get_eos();
  *anyerr = 0;
  for (int i=0; i<npoints; ++i) 
  {
    auto s = state_rho_eps_ye(eos, rho[i], eps[i], ye[i], keyerr[i]);

    if (s) 
    {
      press[i]  = s.press();
      CCTK_REAL cs = s.csnd();
      cs2[i]    = cs*cs;
    }
    else 
    {
      press[i]  = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      cs2[i]    = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      *anyerr   = 1;
    }
  }
}

void eps_cs2_from_rho(const CCTK_INT npoints, 
  const CCTK_REAL* rho, CCTK_REAL* eps, CCTK_REAL* cs2,
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  const eos_barotr& eos = global_eos_cold::get_eos();
  *anyerr = 0;
  for (int i=0; i<npoints; ++i) 
  {
    auto s = state_rho(eos, rho[i], keyerr[i]);
    if (s) 
    {
      eps[i]    = s.eps();
      CCTK_REAL cs = s.csnd();
      cs2[i]    = cs*cs;
    }
    else 
    {
      eps[i]    = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      cs2[i]    = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      *anyerr   = 1;
    }
  }
}



void eps_cs2_from_rho_temp_ye(const CCTK_INT npoints, 
  const CCTK_REAL* rho, CCTK_REAL* eps, const CCTK_REAL* temp, 
  const CCTK_REAL* ye, CCTK_REAL* cs2,
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{  
  const eos_thermal& eos = global_eos_thermal::get_eos();
  *anyerr = 0;
  for (int i=0; i<npoints; ++i) 
  {
    auto s = state_rho_temp_ye(eos, rho[i], temp[i], ye[i], keyerr[i]);
    if (s) 
    {
      eps[i]    = s.eps();
      CCTK_REAL cs = s.csnd();
      cs2[i]    = cs*cs;
    }
    else 
    {
      eps[i]    = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      cs2[i]    = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      *anyerr   = 1;
    }
  }
}


void cs2_from_rho_eps_ye(const CCTK_INT npoints, 
  const CCTK_REAL* rho, const CCTK_REAL* eps, 
  const CCTK_REAL* ye, CCTK_REAL* cs2,
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  const eos_thermal& eos = global_eos_thermal::get_eos();
  *anyerr = 0;
  for (int i=0; i<npoints; ++i) 
  {
    auto s = state_rho_eps_ye(eos, rho[i], eps[i], ye[i], keyerr[i]);

    if (s) 
    {
      CCTK_REAL cs = s.csnd();
      cs2[i]    = cs*cs;
    }
    else 
    {
      cs2[i]    = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      *anyerr   = 1;
    }
  }
}



void eps_DPressByDEps_from_rho_temp_ye(const CCTK_INT npoints, 
  const CCTK_REAL* rho, CCTK_REAL* eps, const CCTK_REAL* temp,
  const CCTK_REAL* ye, CCTK_REAL* dpress_deps, 
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  
  const eos_thermal& eos = global_eos_thermal::get_eos();
  *anyerr = 0;
  for (int i=0; i<npoints; ++i) 
  {
    auto s = state_rho_temp_ye(eos, rho[i], temp[i], ye[i], keyerr[i]);
    if (s) 
    {
      eps[i]         = s.eps();
      dpress_deps[i] = s.dpress_deps();
    }
    else 
    {
      eps[i]          = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      dpress_deps[i]  = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      *anyerr = 1;
    }
  }
}


void DPressByDEps_from_rho_eps_ye(const CCTK_INT npoints, 
  const CCTK_REAL* rho, const CCTK_REAL* eps, const CCTK_REAL* ye, 
  CCTK_REAL* dpress_deps, CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  const eos_thermal& eos = global_eos_thermal::get_eos();
  *anyerr = 0;
  for (int i=0; i<npoints; ++i) 
  {
    auto s = state_rho_eps_ye(eos, rho[i], eps[i], ye[i], keyerr[i]);
    if (s) 
    {
      dpress_deps[i] = s.dpress_deps();
    }
    else 
    {
      dpress_deps[i]  = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      *anyerr = 1;
    }
  }
}



void eps_DPressByDRho_from_rho_temp_ye(const CCTK_INT npoints, 
  const CCTK_REAL* rho, CCTK_REAL* eps, const CCTK_REAL* temp,
  const CCTK_REAL* ye, CCTK_REAL* dpress_drho,
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  const eos_thermal& eos = global_eos_thermal::get_eos();
  *anyerr = 0;
  for (int i=0; i<npoints; ++i) 
  {
    auto s = state_rho_temp_ye(eos, rho[i], temp[i], ye[i], keyerr[i]);
    if (s) 
    {
      eps[i]          = s.eps();
      dpress_drho[i]  = s.dpress_drho();
    }
    else 
    {
      eps[i]          = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      dpress_drho[i]  = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      *anyerr = 1;
    }
  }
}



void DPressByDRho_from_rho_eps_ye(const CCTK_INT npoints, 
  const CCTK_REAL* rho, CCTK_REAL* eps,
  const CCTK_REAL* ye, CCTK_REAL* dpress_drho,
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  const eos_thermal& eos = global_eos_thermal::get_eos();
  *anyerr = 0;
  for (int i=0; i<npoints; ++i) 
  {
    auto s = state_rho_eps_ye(eos, rho[i], eps[i], ye[i], keyerr[i]);
    if (s) 
    {
      dpress_drho[i] = s.dpress_drho();
    }
    else 
    {
      dpress_drho[i]  = std::numeric_limits<CCTK_REAL>::quiet_NaN();
      *anyerr = 1;
    }
  }
}



} // anonymous


extern "C" CCTK_INT RePrimAnd_EOS_Omni_API_GetHandle(const char* name)
{
  const std::string n(name);
  if (n == "RePrimAnd_Evol_EOS")  
  {
    return eos_key_reprimand_thermal;
  }
  if (n == "RePrimAnd_Initial_EOS")  
  {
    return eos_key_reprimand_cold;
  }
  
  CCTK_WARN(0, "Handle requested for EOS not "
                "implemented in RePrimAnd_EOS_Omni_API");
  return 0;
}


extern "C" void RePrimAnd_EOS_Omni_API_press(CCTK_INT eoskey,
  const CCTK_INT havetemp, CCTK_REAL rf_precision,
  const CCTK_INT npoints, const CCTK_REAL* rho,
  CCTK_REAL* eps, CCTK_REAL* temp,
  const CCTK_REAL* ye, CCTK_REAL* press,
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  try {
    if (eoskey == eos_key_reprimand_thermal) 
    {
      if (havetemp == 1) 
      {
        press_eps_from_rho_temp_ye(npoints, 
            rho, eps, temp, ye, press, keyerr, anyerr);      
      }
      else 
      {
        press_from_rho_eps_ye(npoints, 
            rho, eps, ye, press, keyerr, anyerr);    
      }
    }
    else if (eoskey == eos_key_reprimand_cold)
    {
      press_eps_from_rho(npoints, rho, eps, press, keyerr, anyerr);
    } 
    else {
      warn_invalid_handle(eoskey);
    }
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
}

extern "C" void RePrimAnd_EOS_Omni_API_pressOMP(CCTK_INT eoskey,
  const CCTK_INT havetemp, CCTK_REAL rf_precision,
  const CCTK_INT npoints, const CCTK_REAL* rho,
  CCTK_REAL* eps, CCTK_REAL* temp,
  const CCTK_REAL* ye, CCTK_REAL* press,
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  //Not implementing OpenMP version for now
  RePrimAnd_EOS_Omni_API_press(eoskey, havetemp, rf_precision,
      npoints, rho, eps, temp, ye, press, keyerr, anyerr);
}


extern "C" void RePrimAnd_EOS_Omni_API_press_cs2(CCTK_INT eoskey, 
  const CCTK_INT havetemp, CCTK_REAL rf_precision, 
  const CCTK_INT npoints, const CCTK_REAL* rho, 
  CCTK_REAL* eps, CCTK_REAL* temp, const CCTK_REAL* ye,
  CCTK_REAL* press, CCTK_REAL* cs2,
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  try {
    if (eoskey == eos_key_reprimand_thermal) 
    {
      if (havetemp == 1) 
      {
        press_eps_cs2_from_rho_temp_ye(npoints, 
            rho, eps, temp, ye, press, cs2, keyerr, anyerr);
      }
      else 
      {
        press_cs2_from_rho_eps_ye(npoints, 
            rho, eps, ye, press, cs2, keyerr, anyerr);
      }
    }
    else if (eoskey == eos_key_reprimand_cold)
    {
      press_eps_cs2_from_rho(npoints, 
          rho, eps, press, cs2, keyerr, anyerr);
    } 
    else {
      warn_invalid_handle(eoskey);
    }
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
}


extern "C" void RePrimAnd_EOS_Omni_API_cs2(CCTK_INT eoskey,
  const CCTK_INT havetemp, CCTK_REAL rf_precision, 
  const CCTK_INT npoints, const CCTK_REAL* rho,
  CCTK_REAL* eps, CCTK_REAL* temp, const CCTK_REAL* ye,
  CCTK_REAL* cs2, CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  try {
    if (eoskey == eos_key_reprimand_thermal) 
    {
      if (havetemp == 1) 
      {
        eps_cs2_from_rho_temp_ye(npoints, 
            rho, eps, temp, ye, cs2, keyerr, anyerr);
      }
      else 
      {
        cs2_from_rho_eps_ye(npoints, 
            rho, eps, ye, cs2, keyerr, anyerr);
      }
    }
    else if (eoskey == eos_key_reprimand_cold)
    {
      eps_cs2_from_rho(npoints, rho, eps, cs2, keyerr, anyerr);
    } 
    else {
      warn_invalid_handle(eoskey);
    }
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
}


extern "C" void RePrimAnd_EOS_Omni_API_DPressByDEps(CCTK_INT eoskey,
  const CCTK_INT havetemp, CCTK_REAL rf_precision, 
  const CCTK_INT npoints, const CCTK_REAL* rho, 
  CCTK_REAL* eps, CCTK_REAL* temp,
  const CCTK_REAL* ye, CCTK_REAL* dpress_deps, 
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  try {
    if (eoskey == eos_key_reprimand_thermal) 
    {
      if (havetemp == 1) 
      {
        eps_DPressByDEps_from_rho_temp_ye(npoints, 
            rho, eps, temp, ye, dpress_deps, keyerr, anyerr);
        
      }
      else 
      {
        DPressByDEps_from_rho_eps_ye(npoints, 
            rho, eps, ye, dpress_deps, keyerr, anyerr);
      }
    }
    else if (eoskey == eos_key_reprimand_cold)
    {
      CCTK_WARN(0, "Partial derivative of pressure w.r.t epsilon "
                    "requested for a barotropic EOS");
    } 
    else {
      warn_invalid_handle(eoskey);
    }
      
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
}

extern "C" void RePrimAnd_EOS_Omni_API_DPressByDRho(CCTK_INT eoskey,
  const CCTK_INT havetemp, CCTK_REAL rf_precision, 
  const CCTK_INT npoints, const CCTK_REAL* rho, 
  CCTK_REAL* eps, CCTK_REAL* temp,
  const CCTK_REAL* ye, CCTK_REAL* dpress_drho,
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  try {
    if (eoskey == eos_key_reprimand_thermal) 
    {
      if (havetemp == 1) 
      {
        eps_DPressByDRho_from_rho_temp_ye(npoints, 
            rho, eps, temp, ye, dpress_drho, keyerr, anyerr);
      }
      else 
      {
        DPressByDRho_from_rho_eps_ye(npoints, 
            rho, eps, ye, dpress_drho, keyerr, anyerr);
      }
    }
    else if (eoskey == eos_key_reprimand_cold)
    {
      CCTK_WARN(0, "Partial derivative of pressure w.r.t rho "
                    "requested for a barotropic EOS");
    } 
    else {
      warn_invalid_handle(eoskey);
    }
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
}




extern "C" void RePrimAnd_EOS_Omni_API_dpderho_dpdrhoe(
  CCTK_INT eoskey, const CCTK_INT havetemp, CCTK_REAL rf_precision,
  const CCTK_INT npoints, const CCTK_REAL* rho,
  CCTK_REAL* eps, CCTK_REAL* temp, const CCTK_REAL* ye,
  CCTK_REAL* dpderho, CCTK_REAL* dpdrhoe,
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  CCTK_WARN(0, "EOS_Omni_dpderho_dpdrhoe not implemented by "
               "RePrimAnd_EOS_Omni_API");
}

extern "C" void RePrimAnd_EOS_Omni_API_EpsFromPress(CCTK_INT eoskey,
  const CCTK_INT havetemp, CCTK_REAL rf_precision,
  const CCTK_INT npoints, const CCTK_REAL* rho,
  CCTK_REAL* eps, CCTK_REAL* temp, const CCTK_REAL* ye,
  const CCTK_REAL* press, CCTK_REAL* xeps, 
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  CCTK_WARN(0, "EOS_Omni_EpsFromPress not implemented by " 
               "RePrimAnd_EOS_Omni_API");
}


extern "C" void RePrimAnd_EOS_Omni_API_RhoFromPressEpsTempEnt(
  CCTK_INT eoskey, const CCTK_INT havetemp, CCTK_REAL rf_precision,
  const CCTK_INT npoints, CCTK_REAL* rho, CCTK_REAL* eps,
  CCTK_REAL* temp, CCTK_REAL* ent, const CCTK_REAL* ye, 
  const CCTK_REAL* press, CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  CCTK_WARN(0, "EOS_Omni_RhoFromPressEpsTempEnt not implemented by "
               "RePrimAnd_EOS_Omni_API");
}


extern "C" void RePrimAnd_EOS_Omni_API_PressEpsTempYe_from_Rho(
  CCTK_INT eoskey, const CCTK_INT havetemp, CCTK_REAL rf_precision,
  const CCTK_INT npoints, const CCTK_REAL* rho, CCTK_REAL* eps,
  CCTK_REAL* temp, CCTK_REAL* ye, CCTK_REAL* press, 
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  CCTK_WARN(0, "EOS_Omni_PressEpsTempYe_from_Rho not implemented by "
               "RePrimAnd_EOS_Omni_API");
}


extern "C" void RePrimAnd_EOS_Omni_API_short(CCTK_INT eoskey,  
  const CCTK_INT havetemp, CCTK_REAL rf_precision, 
  const CCTK_INT npoints, const CCTK_REAL* rho, 
  CCTK_REAL* eps, CCTK_REAL* temp, const CCTK_REAL* ye,
  CCTK_REAL* press, CCTK_REAL* entropy, CCTK_REAL* cs2,
  CCTK_REAL* dedt, CCTK_REAL* dpderho, CCTK_REAL* dpdrhoe,
  CCTK_REAL* munu, CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  CCTK_WARN(0, "EOS_Omni_short not implemented by " 
               "RePrimAnd_EOS_Omni_API");
}



extern "C" void RePrimAnd_EOS_Omni_API_full(CCTK_INT eoskey, 
  const CCTK_INT havetemp, CCTK_REAL rf_precision,
  const CCTK_INT npoints, const CCTK_REAL* rho,
  CCTK_REAL* eps, CCTK_REAL* temp, const CCTK_REAL* ye,
  CCTK_REAL* press, CCTK_REAL* entropy, CCTK_REAL* cs2,
  CCTK_REAL* dedt, CCTK_REAL* dpderho, CCTK_REAL* dpdrhoe, 
  CCTK_REAL* xa, CCTK_REAL* xh, CCTK_REAL* xn, CCTK_REAL* xp, 
  CCTK_REAL* abar, CCTK_REAL* zbar, CCTK_REAL* mue, 
  CCTK_REAL* mun, CCTK_REAL* mup, CCTK_REAL* muhat,
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  CCTK_WARN(0, "EOS_Omni_full not implemented by "
               "RePrimAnd_EOS_Omni_API");
}


extern "C" void RePrimAnd_EOS_Omni_API_DEpsByDRho_DEpsByDPress(
  CCTK_INT eoskey, const CCTK_INT havetemp, CCTK_REAL rf_precision,
  const CCTK_INT npoints, const CCTK_REAL* rho, CCTK_REAL* eps,
  CCTK_REAL* temp, const CCTK_REAL* ye, CCTK_REAL* DEpsByDRho, 
  CCTK_REAL* DEpsByDPress, CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  CCTK_WARN(0, "EOS_Omni_DEpsByDRho_DEpsByDPress not implemented "
               "by RePrimAnd_EOS_Omni_API");
}



extern "C" void RePrimAnd_EOS_Omni_API_press_f_hrho_v2_rhoW(
  CCTK_INT eoskey, const CCTK_INT havetemp, CCTK_REAL rf_precision,
  const CCTK_INT npoints, const CCTK_REAL* hrho,
  const CCTK_REAL* v2, const CCTK_REAL* rhoW, CCTK_REAL* eps,
  CCTK_REAL* temp, const CCTK_REAL* ye, CCTK_REAL* press,
  CCTK_INT* keyerr,  CCTK_INT* anyerr)
{
  CCTK_WARN(0, "EOS_Omni_press_f_hrho_v2_rhoW not implemented by "
               "RePrimAnd_EOS_Omni_API");
}




extern "C" void RePrimAnd_EOS_Omni_API_dpdhrho_f_hrho_v2_rhoW(
  CCTK_INT eoskey, const CCTK_INT havetemp, CCTK_REAL rf_precision,
  const CCTK_INT npoints, const CCTK_REAL* hrho, const CCTK_REAL* v2, 
  const CCTK_REAL* rhoW, CCTK_REAL* eps, CCTK_REAL* temp,
  const CCTK_REAL* ye, CCTK_REAL* dpdhrho,
  CCTK_INT* keyerr, CCTK_INT* anyerr)
{
  CCTK_WARN(0, "EOS_Omni_dpdhrho_f_hrho_v2_rhoW not implemented by "
               "RePrimAnd_EOS_Omni_API");
}


extern "C" void RePrimAnd_EOS_Omni_API_dpdv2_f_hrho_v2_rhoW(
  CCTK_INT eoskey, const CCTK_INT havetemp, CCTK_REAL rf_precision,
  const CCTK_INT npoints, const CCTK_REAL* hrho, const CCTK_REAL* v2, 
  const CCTK_REAL* rhoW, CCTK_REAL* eps, CCTK_REAL* temp,
  const CCTK_REAL* ye, CCTK_REAL* dpdv2, 
  CCTK_INT* keyerr, CCTK_INT* anyerr) 
{
  CCTK_WARN(0, "EOS_Omni_dpdv2_f_hrho_v2_rhoW not implemented by "
               "RePrimAnd_EOS_Omni_API");
}
  
