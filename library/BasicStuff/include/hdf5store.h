#ifndef HDF5STORE_H
#define HDF5STORE_H

#include <string>
#include "datastore.h"

namespace EOS_Toolkit {

datasource make_hdf5_file_source(std::string path);
datasink make_hdf5_file_sink(std::string path);


} // namespace EOS_Toolkit

#endif
