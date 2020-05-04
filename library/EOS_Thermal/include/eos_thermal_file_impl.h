#ifndef EOS_THERMAL_FILE_IMPL_H
#define EOS_THERMAL_FILE_IMPL_H

#include <string>
#include "global_registry.h"
#include "hdf5imple.h"
#include "eos_thermal.h"
#include "unitconv.h"
 

namespace EOS_Toolkit {


namespace implementations {

struct reader_eos_thermal {
  virtual eos_thermal load(const h5grp& f, const units& u) const = 0;
  virtual ~reader_eos_thermal() {} 
};

using registry_reader_eos_thermal = global_registry<reader_eos_thermal>;

} // namespace implementations


} // namespace EOS_Toolkit

#endif

