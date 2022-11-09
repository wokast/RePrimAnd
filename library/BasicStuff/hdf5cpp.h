#ifndef HDF5CPP_H
#define HDF5CPP_H

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
class h5_resource {
  struct raw {
    const hid_t h;
    raw(hid_t h_) : h{h_} {}
    ~raw() {if (h>=0) C::close(h);}
  };

  std::shared_ptr<const raw> p;

  public:
  
  template<class... A>
  explicit h5_resource(const A&... a) 
  : p(std::make_shared<raw>(C::open(a...))) {}
  
  h5_resource(const h5_resource&) = default;
  h5_resource(h5_resource&&) = default;
  h5_resource& operator=(h5_resource&&) = default;
  h5_resource& operator=(const h5_resource&) = default;

  hid_t use() const {
    if (p->h < 0) throw std::runtime_error(C::err_msg());
    return p->h;
  }
  
  explicit operator bool() const {return p->h >= 0;}
};




struct h5api_file_read {
  static hid_t open(std::string path);
  static void close(hid_t h);
  static const char* err_msg();
};

using h5_file_read = h5_resource<h5api_file_read>;

struct h5api_file_write {
  static hid_t open(std::string path);
  static void close(hid_t h);
  static const char* err_msg();
};

using h5_file_write = h5_resource<h5api_file_write>;


struct h5api_group_read {
  static hid_t open(hid_t loc, std::string name) ;
  static hid_t open(const h5_file_read& f, std::string name);
  static hid_t open(const h5_resource<h5api_group_read>& g, 
                    std::string name);
  static void close(hid_t h);
  static const char* err_msg();
};

using h5_group_read = h5_resource<h5api_group_read>;


struct h5api_group_write {
  static hid_t open(hid_t loc, std::string name) ;
  static hid_t open(const h5_file_write& f, std::string name);
  static hid_t open(const h5_resource<h5api_group_write>& g, 
                    std::string name);
  static void close(hid_t h);
  static const char* err_msg();
};

using h5_group_write = h5_resource<h5api_group_write>;




struct h5api_dspc_scalar {
  static hid_t open();
  static void close(hid_t h);
  static const char* err_msg();
};

using h5_dspc_scalar = h5_resource<h5api_dspc_scalar>;


struct h5api_attr_read {
  static hid_t open(const h5_file_read& f, std::string name);
  static hid_t open(const h5_group_read& g, std::string name);
  static void close(hid_t h);
  static const char* err_msg();
};

using h5_attr_read = h5_resource<h5api_attr_read>;



struct h5api_attr_write {
  static hid_t open(hid_t loc, std::string n, 
                    const h5_dspc_scalar& s, hid_t t);
                    
  template<class R, class T>
  static hid_t open(const R& f, std::string n, 
                    const h5_dspc_scalar& s, T&t)
  {
    return open(f.use(), n, s, t);
  }
  
  
  static void close(hid_t h);
  static const char* err_msg();
};

using h5_attr_write = h5_resource<h5api_attr_write>;


struct h5api_dtyp_read {
  static hid_t open(const h5_attr_read& a);
  static void close(hid_t h);
  static const char* err_msg();
};

using h5_dtyp_read = h5_resource<h5api_dtyp_read>;


struct h5api_dtyp_string {
  static hid_t open();
  static void close(hid_t h);
  static const char* err_msg();
};

using h5_dtyp_string = h5_resource<h5api_dtyp_string>;




struct h5api_dset_read {
  static hid_t open(hid_t loc, std::string name);
  static hid_t open(const h5_file_read& f, std::string name);
  static hid_t open(h5_group_read g, std::string name);
  static void close(hid_t h);
  static const char* err_msg();
};

using h5_dset_read = h5_resource<h5api_dset_read>;



struct h5api_dspc_write {
  template<std::size_t N>
  static hid_t open(const std::array<hsize_t, N> e)
  {
    return H5Screate_simple(N, &(e[0]), NULL);
  }
  static void close(hid_t h);
  static const char* err_msg();
};

using h5_dspc_write = h5_resource<h5api_dspc_write>;



struct h5api_dset_write {
  static hid_t open(hid_t loc, std::string n, hid_t t, 
                    const h5_dspc_write& s);
  template<class R>
  static hid_t open(const R& f, std::string n,  hid_t t,
                    const h5_dspc_write& s)
  {
    return open(f.use(), n, t, s);
  }
  
  static void close(hid_t h);
  static const char* err_msg();
};

using h5_dset_write = h5_resource<h5api_dset_write>;


struct h5api_dspc_read {
  static hid_t open(const h5_dset_read &ds);
  static void close(hid_t h);
  static const char* err_msg();
};

using h5_dspc_read = h5_resource<h5api_dspc_read>;






template<class R>
bool has_link(const R&r, std::string n)
{
  return (H5Lexists(r.use(), n.c_str(), H5P_DEFAULT) > 0);
}

template<class R>
bool has_attrib(const R&r, std::string n)
{
  return (H5Aexists(r.use(), n.c_str()) > 0);
}

template<class R>
bool has_entry(const R&r, std::string n)
{
  return (has_link(r, n) || has_attrib(r,n));
}

template<class T>
void read_attr(const h5_attr_read& a, T& d)
{
  if (H5Aread(a.use(), detail::h5_types<T>::id(), &d) < 0) {
    throw std::runtime_error("HDF5: problem reading attribute");
  }
}



template<class R, class T>
void read_attr(const R& r, std::string n, T& d)
{
  h5_attr_read a(r, n);
  read_attr(a, d);
}


template<class R>
void read_attr(const R& r, std::string n, bool& d)
{
  int i;
  read_attr(r, n, i);
  d = (i != 0);
}



template<class T>
void write_attr(const h5_attr_write& a, hid_t t, T& d)
{
  if (H5Awrite(a.use(), t, &d) < 0) {
    throw std::runtime_error("HDF5: problem writing attribute");
  }
}

template<class R, class T>
void write_attr(R& r, std::string n, hid_t t, const T& v)
{
  h5_dspc_scalar s;
  h5_attr_write a(r, n, s, t);
  write_attr(a, t, v);
}


template<class R, class T>
void write_attr(R& r, std::string n, const T& v)
{
  hid_t t = detail::h5_types<T>::id();
  write_attr(r, n, t, v);
}

template<class R>
void write_attr(R& r, std::string n, const std::string& v)
{
  h5_dtyp_string t;
  write_attr(r, n, t.use(), v);
}


template<class R>
void write_attr(R& r, std::string n, const bool& v)
{
  int i = v ? -1 : 0;
  write_attr(r, n, i);
}



template<std::size_t N>
std::array<hsize_t, N> get_extent(const h5_dspc_read& dspc) 
{
  std::array<hsize_t, N> ext;
  if (N != H5Sget_simple_extent_ndims(dspc.use())) 
  {
    throw std::runtime_error(
               "HDF5: dataset with unexpected number dimensions.");
  }
  if (N != H5Sget_simple_extent_dims(dspc.use(), &(ext[0]), nullptr)) 
  {
    throw std::runtime_error(
               "HDF5: problem getting dataset extent.");
  }
  return ext;
}

template<std::size_t N>
hsize_t extent2size(const std::array<hsize_t, N>& ext)
{
  hsize_t sz = 1;
  for (auto s : ext) sz *= s;
  return sz;
}



template<class T>
void read_data(const h5_dset_read &dset, T* buf, std::size_t size)
{
  h5_dspc_read dspc(dset);
  hssize_t npts = H5Sget_simple_extent_npoints(dspc.use());
  if (0 > npts) {
    throw std::runtime_error("HDF5: problem getting data size");
  }
  if (std::size_t(npts) != size) {
    throw std::runtime_error("HDF5: unexpected dataset size");
  }
  
  
  if (0 > H5Dread(dset.use(), detail::h5_types<T>::id(), 
                  H5S_ALL, dspc.use(), H5P_DEFAULT, buf) )
  {
    throw std::runtime_error("HDF5: problem reading dataset");
  }    
}

template<class T>
void read_data(const h5_dset_read& dset, std::vector<T>& v)
{
  h5_dspc_read dspc(dset);
  auto ext = get_extent<1>(dspc);
  v.resize(extent2size(ext));
  read_data(dset, v.data(), v.size());
}

template<class R, class T>
void read_data(const R& r, std::string n, std::vector<T>& v)
{
  h5_dset_read dset(r, n);
  read_data(dset, v);
}




template<class T>
void write_data(const h5_dset_write& dset, const std::vector<T>& v)
{
  if (H5Dwrite(dset.use(), detail::h5_types<T>::id(), 
           H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) < 0 ) 
  {
    throw std::runtime_error("HDF5: problem writing dataset");
  }
}

template<class R, class T>
void write_data(const R& r, std::string n, const std::vector<T>& v)
{
  std::array<hsize_t, 1> ext{{v.size()}}; 
  h5_dspc_write s(ext);
  h5_dset_write d(r, n, detail::h5_types<T>::id(), s);
  write_data(d, v);
}


void read_attr(const h5_attr_read& a, std::string& d);

} //namespace detail



} // namespace EOS_Toolkit

#endif
