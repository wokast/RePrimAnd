#ifndef REPRIMAND_GLOBAL_EOS_H
#define REPRIMAND_GLOBAL_EOS_H

#include "reprimand/eos_thermal.h"
#include "reprimand/eos_barotropic.h"


namespace RePrimAnd_Global_EOS {

class global_eos_thermal {
  static EOS_Toolkit::eos_thermal eos;
  static bool initialized;
  public:
  static const EOS_Toolkit::eos_thermal& get_eos();
  static void set_eos(const EOS_Toolkit::eos_thermal& eos_);
};

class global_eos_cold {
  static EOS_Toolkit::eos_barotr eos;
  static bool initialized;
  public:
  static const EOS_Toolkit::eos_barotr& get_eos();
  static void set_eos(const EOS_Toolkit::eos_barotr& eos_);
};

}

#endif
