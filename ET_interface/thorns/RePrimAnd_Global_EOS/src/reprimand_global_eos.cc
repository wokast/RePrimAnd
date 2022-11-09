#include "reprimand_global_eos.h"
#include <stdexcept>

using namespace EOS_Toolkit;

namespace RePrimAnd_Global_EOS {
  
eos_thermal global_eos_thermal::eos;
bool        global_eos_thermal::initialized = false;

eos_barotr  global_eos_cold::eos;
bool        global_eos_cold::initialized = false;

const eos_thermal& global_eos_thermal::get_eos()
{
  if (!initialized) {
    throw std::logic_error("global_eos_thermal: access to eos requested before initialization");
  }
  return eos;
}

void global_eos_thermal::set_eos(const eos_thermal& eos_)
{
  if (initialized) {
    throw std::logic_error("global_eos_thermal: multiple initialization attempted");
  }
  eos         = eos_;
  initialized = true;
}

const eos_barotr& global_eos_cold::get_eos()
{
  if (!initialized) {
    throw std::logic_error("global_eos_cold: access to eos requested before initialization");
  }
  return eos;
}

void global_eos_cold::set_eos(const eos_barotr& eos_)
{
  if (initialized) {
    throw std::logic_error("global_eos_cold: multiple initialization attempted");
  }
  eos         = eos_;
  initialized = true;
}

}

