#include "hdf5imple.h"

namespace EOS_Toolkit {


void read(const h5attr& a, bool& b)
{
  int i;
  read(a,i);
  b = (i != 0);
}

void read(const h5attr& a, std::string& d)
{
  detail::h5dtype t(a);
  if (H5Tget_class(t.use()) != H5T_STRING) {
    throw std::runtime_error("HDF5: expected string attribute");
  }

  if (H5Tis_variable_str(t.use()) <= 0) {
    throw std::runtime_error("HDF5: expected variable length string");
  }
  
  char* buf = nullptr;
  
  if (H5Aread(a.use(), t.use(), &buf) < 0) {
    throw std::runtime_error("HDF5: problem reading attribute");
  }
  
  d = std::string(buf);
  
  H5free_memory(buf);
}


}
