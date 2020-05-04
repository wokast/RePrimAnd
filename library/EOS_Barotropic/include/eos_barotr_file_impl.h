#ifndef EOS_BAROTR_FILE_IMPL_H
#define EOS_BAROTR_FILE_IMPL_H

#include <string>
#include "global_registry.h"
#include "hdf5imple.h"
#include "eos_barotropic.h"
#include "unitconv.h"
 

namespace EOS_Toolkit {


namespace implementations {

struct reader_eos_barotr {
  virtual eos_barotr load(const h5grp& f, const units& u) const = 0;
  virtual ~reader_eos_barotr() {} 
};

using registry_reader_eos_barotr = global_registry<reader_eos_barotr>;

} // namespace implementations


} // namespace EOS_Toolkit

#endif

