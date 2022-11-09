#include "hdf5cpp.h"
#include <cassert>



namespace EOS_Toolkit {

namespace detail {

void read_attr(const h5_attr_read& a, std::string& d)
{
  h5_dtyp_read t(a);
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
  assert(buf);
  
  d = std::string(buf);
  
  H5free_memory(buf);
}



hid_t h5api_file_read::open(std::string path) 
{
  return H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
}

void h5api_file_read::close(hid_t h) 
{
  H5Fclose(h);
}

const char* h5api_file_read::err_msg() 
{
  return "HDF5: problem reading file";
}



hid_t h5api_file_write::open(std::string path) 
{
  return H5Fcreate(path.c_str(), H5F_ACC_EXCL, 
                   H5P_DEFAULT, H5P_DEFAULT);
}


void h5api_file_write::close(hid_t h) 
{
  H5Fclose(h);
}

const char* h5api_file_write::err_msg() 
{
  return "HDF5: problem creating file";
}




hid_t h5api_attr_read::open(const h5_file_read& f, std::string name) 
{
  return H5Aopen(f.use(), name.c_str(), H5P_DEFAULT);
}

hid_t h5api_attr_read::open(const h5_group_read& g, std::string name) 
{
  return H5Aopen(g.use(), name.c_str(), H5P_DEFAULT);
}

void h5api_attr_read::close(hid_t h) 
{
  H5Aclose(h);
}

const char* h5api_attr_read::err_msg() 
{
  return "HDF5: problem finding attribute";
}



hid_t h5api_attr_write::open(hid_t loc, std::string n, 
                    const h5_dspc_scalar& s, hid_t t) 
{
  return H5Acreate2(loc, n.c_str(), t, s.use(),  
                    H5P_DEFAULT, H5P_DEFAULT);
}

void h5api_attr_write::close(hid_t h) 
{
  H5Aclose(h);
}

const char* h5api_attr_write::err_msg() 
{
  return "HDF5: problem creating attribute";
}



hid_t h5api_dtyp_read::open(const h5_attr_read& a) 
{
  return H5Aget_type(a.use());
}

void h5api_dtyp_read::close(hid_t h) 
{
  H5Tclose(h);
}

const char* h5api_dtyp_read::err_msg() 
{
  return "HDF5: problem obtaining attribute data type";
}


hid_t h5api_dset_read::open(hid_t loc, std::string name) 
{
  if (H5Lexists(loc, name.c_str(), H5P_DEFAULT) <= 0) return -1;
  return H5Dopen(loc, name.c_str(), H5P_DEFAULT);
}

hid_t h5api_dset_read::open(const h5_file_read& f, std::string name) 
{
  return open(f.use(), name);
}

hid_t h5api_dset_read::open(h5_group_read g, std::string name) 
{
  return open(g.use(), name);
}

void h5api_dset_read::close(hid_t h) 
{
  H5Dclose(h);
}

const char* h5api_dset_read::err_msg() 
{
  return "HDF5: problem opening data set";
}


hid_t h5api_dset_write::open(hid_t loc, std::string n,  hid_t t, 
                             const h5_dspc_write& s) 
{
  return H5Dcreate2(loc, n.c_str(), t, s.use(), 
                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

void h5api_dset_write::close(hid_t h) 
{
  H5Dclose(h);
}

const char* h5api_dset_write::err_msg() 
{
  return "HDF5: problem creating dataset";
}



hid_t h5api_dspc_read::open(const h5_dset_read &ds) 
{ 
  return H5Dget_space(ds.use());
}

void h5api_dspc_read::close(hid_t h) 
{
  H5Sclose(h);
}

const char* h5api_dspc_read::err_msg() 
{
  return "HDF5: problem obtaining dataset data space";
}



void h5api_dspc_write::close(hid_t h) 
{
  H5Sclose(h);
}

const char* h5api_dspc_write::err_msg() 
{
  return "HDF5: problem creating data space";
}



hid_t h5api_dspc_scalar::open() 
{ 
  return H5Screate(H5S_SCALAR);
}

void h5api_dspc_scalar::close(hid_t h) 
{
  H5Sclose(h);
}

const char* h5api_dspc_scalar::err_msg() 
{
  return "HDF5: problem creating scalar data space";
}


hid_t h5api_dtyp_string::open() 
{ 
  hid_t t = H5Tcopy(H5T_C_S1);
  if (H5Tset_size(t, H5T_VARIABLE) < 0) return -1;
  return t;

}

void h5api_dtyp_string::close(hid_t h) 
{
  //~ H5Sclose(h);
}

const char* h5api_dtyp_string::err_msg() 
{
  return "HDF5: problem creating string data type";
}




hid_t h5api_group_read::open(hid_t loc, std::string name) 
{
  if (H5Lexists(loc, name.c_str(), H5P_DEFAULT) <= 0) return -1;
  return H5Gopen(loc, name.c_str(), H5P_DEFAULT);
}
  
hid_t h5api_group_read::open(const h5_file_read& f, std::string name) 
{
  return open(f.use(), name);
}
  
hid_t h5api_group_read::open(const h5_resource<h5api_group_read>& g, 
                    std::string name) 
{
  return open(g.use(), name);
}

void h5api_group_read::close(hid_t h) 
{
  H5Gclose(h);
}

const char* h5api_group_read::err_msg() 
{
  return "HDF5: problem opening group";
}





hid_t h5api_group_write::open(hid_t loc, std::string name) 
{
  return H5Gcreate(loc, name.c_str(), 
                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}
  
hid_t h5api_group_write::open(const h5_file_write& f, std::string name) 
{
  return open(f.use(), name);
}
  
hid_t h5api_group_write::open(const h5_resource<h5api_group_write>& g, 
                    std::string name) 
{
  return open(g.use(), name);
}

void h5api_group_write::close(hid_t h) 
{
  H5Gclose(h);
}

const char* h5api_group_write::err_msg() 
{
  return "HDF5: problem creating group";
}




} // namespace detail 

} // namespace EOS_Toolkit 
