#ifndef HDF5IMPLE_H
#define HDF5IMPLE_H

#include <hdf5.h>
#include <memory>
#include <string>
#include <vector>
#include <array>
#include <stdexcept>

namespace EOS_Toolkit {


namespace detail {

template<class T> struct h5_types;
template<> struct h5_types<float> {
  static hid_t id() {return H5T_NATIVE_FLOAT;}
};
template<> struct h5_types<double> {
  static hid_t id() {return H5T_NATIVE_DOUBLE;}
};
template<> struct h5_types<int> {
  static hid_t id() {return H5T_NATIVE_INT;}
};
template<> struct h5_types<unsigned int> {
  static hid_t id() {return H5T_NATIVE_UINT;}
};


template<class C>
class h5id {
  struct raw {
    const hid_t h;
    raw(hid_t h_) : h{h_} {}
    ~raw() {if (h>=0) C::close(h);}
  };

  std::shared_ptr<const raw> p;

  public:
  
  template<class... A>
  explicit h5id(const A&... a) 
  : p(std::make_shared<raw>(C::open(a...))) {}
  
  h5id(const h5id&) = default;
  h5id(h5id&&) = default;
  h5id& operator=(h5id&&) = default;
  h5id& operator=(const h5id&) = default;

  hid_t use() const {
    if (p->h < 0) throw std::runtime_error(C::err_msg());
    return p->h;
  }
  
  explicit operator bool() const {return p->h >= 0;}
};


struct h5id_file {
  static hid_t open(std::string path) 
  {
    return H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  }
  static void close(hid_t h) {H5Fclose(h);}
  static const char* err_msg() {return "HDF5: invalid file";}
};

using h5file = h5id<h5id_file>;

class h5id_grp {
  static hid_t open(hid_t loc, std::string name) 
  {
    if (H5Lexists(loc, name.c_str(), H5P_DEFAULT) <= 0) return -1;
    return H5Gopen(loc, name.c_str(), H5P_DEFAULT);
  }
  public:
  static hid_t open(const h5file& f, std::string name) 
  {
    return open(f.use(), name);
  }
  static hid_t open(const h5id<h5id_grp>& g, std::string name) 
  {
    return open(g.use(), name);
  }
  static void close(hid_t h) {H5Gclose(h);}
  static const char* err_msg() {return "HDF5: invalid group";}
};

using h5grp = h5id<h5id_grp>;



class h5id_dset {
  static hid_t open(hid_t loc, std::string name) 
  {
    if (H5Lexists(loc, name.c_str(), H5P_DEFAULT) <= 0) return -1;
    return H5Dopen(loc, name.c_str(), H5P_DEFAULT);
  }
  public:
  static hid_t open(h5file f, std::string name) {
    return open(f.use(), name);
  }
  static hid_t open(h5grp g, std::string name) {
    return open(g.use(), name);
  }
  static void close(hid_t h) {H5Dclose(h);}
  static const char* err_msg() {return "HDF5: invalid dataset";}
};

using h5dset = h5id<h5id_dset>;


struct h5id_dspc {
  static hid_t open(const h5dset &ds) { 
    return H5Dget_space(ds.use());
  }
  static void close(hid_t h) {H5Sclose(h);}
  static const char* err_msg() {return "HDF5: invalid data space";}
};

using h5dspc = h5id<h5id_dspc>;


struct h5id_attr {
  static hid_t open(const h5file& f, std::string name) 
  {
    return H5Aopen(f.use(), name.c_str(), H5P_DEFAULT);
  }
  static hid_t open(const h5grp& g, std::string name) 
  {
    return H5Aopen(g.use(), name.c_str(), H5P_DEFAULT);
  }
  static hid_t open(const h5dset& d, std::string name) 
  {
    return H5Aopen(d.use(), name.c_str(), H5P_DEFAULT);
  }
  static void close(hid_t h) {H5Aclose(h);}
  static const char* err_msg() {return "HDF5: invalid attribute";}
};

using h5attr = h5id<h5id_attr>;

struct h5id_dtype {
  static hid_t open(const h5attr& a) 
  {
    return H5Aget_type(a.use());
  }
  static hid_t open(const h5dset& d) 
  {
    return H5Dget_type(d.use());
  }
  static void close(hid_t h) {H5Tclose(h);}
  static const char* err_msg() {return "HDF5: invalid data type";}
};

using h5dtype = h5id<h5id_dtype>;


}


using detail::h5file;
using detail::h5grp;
using detail::h5dset;
using detail::h5dspc;
using detail::h5attr;


template<class T>
class h5buf {
  T* p;
  hsize_t sz;
  
  public:
  
  h5buf()             = delete;
  h5buf(const h5buf&) = default;
  h5buf(h5buf&&)      = default;
  h5buf& operator=(const h5buf&) = default;
  h5buf& operator=(h5buf&&)      = default;

  h5buf(T* p_, hsize_t sz_) : p{p_}, sz{sz_}
  {
    if (p == nullptr)
      throw std::logic_error("HDF5: invalid memory destination");
  }
  h5buf(std::vector<T>& v) : h5buf(v.data(), v.size()) {}
  
  T* get(hsize_t minsz) const 
  {
    if (minsz > sz) 
      throw std::runtime_error("HDF5: insufficient buffer size");
    return p;
  }
};

template<hsize_t N>
class h5ext {
  using ext_t = std::array<hsize_t,N>;
  
  ext_t ext; 
  hsize_t tsz;
  
  public:
  h5ext(const h5ext&) = default;
  h5ext(h5ext&&) = default;
  h5ext& operator=(const h5ext&) = default;
  h5ext& operator=(h5ext&&) = default;

  h5ext(const h5dspc& dspc) {
    if (N != H5Sget_simple_extent_ndims(dspc.use())) 
    {
      throw std::runtime_error(
                 "HDF5: wrong number of dataset dimensions.");
    }
    if (N != H5Sget_simple_extent_dims(dspc.use(), &(ext[0]), nullptr)) 
    {
      throw std::runtime_error(
                 "HDF5: problem getting dataset extent.");
    }
    tsz = 1;
    for (auto s : ext) tsz *= s;
  }
  
  h5ext(const h5dset& ds) : h5ext{h5dspc(ds)} {}
  

  ext_t extent() const {return ext;}
  hsize_t size() const {return tsz;}
};


template<class T>
void read(const h5dset &dset, const h5buf<T>& buf)
{
  h5dspc dspc(dset);
  hssize_t npts = H5Sget_simple_extent_npoints(dspc.use());
  if (npts<0) {
    throw std::runtime_error("HDF5: problem getting data size");
  }
  
  if (0 > H5Dread(dset.use(), detail::h5_types<T>::id(), 
                  H5S_ALL, dspc.use(), H5P_DEFAULT, buf.get(npts)) )
  {
    throw std::runtime_error("HDF5: problem reading dataset");
  }    
}

template<class T>
void read(const h5dset &dset, std::vector<T>& v) 
{
  h5ext<1> ext{dset};
  v.resize(ext.size());
  h5buf<T> b(v);
  read(dset, b);
}
  


template<class T>
void read(const h5attr& a, T& d)
{
  if (H5Aread(a.use(), detail::h5_types<T>::id(), &d) < 0) {
    throw std::runtime_error("HDF5: problem reading attribute");
  }
}

void read(const h5attr& a, bool& b);

void read(const h5attr& a, std::string& d);


template<class T>
std::vector<T> get_vector(const h5dset& d)
{
  std::vector<T> v;
  read(d, v);
  return v;
}

template<class T>
std::vector<T> get_vector(const h5grp& g, std::string name)
{
  return get_vector<T>(h5dset(g, name));
}

template<class T>
std::vector<T> get_vector(const h5file& f, std::string name)
{
  
  return get_vector<T>(h5dset(f, name));
}

template<class T>
T get_attribute(const h5attr& a)
{
  T t;
  read(a, t);
  return t;
}

template<class T>
T get_attribute(h5file f, std::string name)
{
  return get_attribute<T>(h5attr(f, name));
}

template<class T>
T get_attribute(h5grp g, std::string name)
{
  return get_attribute<T>(h5attr(g, name));
}

 

} // namespace EOS_Toolkit

#endif
