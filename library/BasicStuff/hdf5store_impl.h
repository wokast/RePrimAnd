#ifndef HDF5STORE_IMPL_H
#define HDF5STORE_IMPL_H

//~ #include <hdf5.h>
#include <memory>
#include <string>
#include <vector>
#include <array>

#include "datastore.h"
#include "hdf5cpp.h"

namespace EOS_Toolkit {


namespace detail {


class h5group_source_impl : public source_impl {
  const h5_group_read res;
  
  public:
  h5group_source_impl(const h5_group_read res_);
  
  bool has_data(std::string) const final;
  void read(std::string, double&) const final;
  void read(std::string, int&) const final;
  void read(std::string, bool&) const final;
  void read(std::string, std::string&) const final;
  void read(std::string, std::vector<double>&) const final;
  void read(std::string, std::vector<int>&) const final;
  
  bool has_group(std::string) const final;
  std::shared_ptr<source_impl> group(std::string) const final;
  
  ~h5group_source_impl() = default;
};


class h5file_source_impl : public source_impl {
  const h5_file_read res;
  
  public:
  h5file_source_impl(std::string path);
  bool has_data(std::string) const final;
  void read(std::string, double&) const final;
  void read(std::string, int&) const final;
  void read(std::string, bool&) const final;
  void read(std::string, std::string&) const final;
  void read(std::string, std::vector<double>&) const final;
  void read(std::string, std::vector<int>&) const final;
  
  bool has_group(std::string) const final;
  std::shared_ptr<source_impl> group(std::string) const final;
  
  ~h5file_source_impl() = default;
};


class h5file_sink_impl : public sink_impl {
  const h5_file_write res;
  
  public:
  h5file_sink_impl(std::string path);
  void write(std::string, const double&) final;
  void write(std::string, const int&) final;
  void write(std::string, const bool&) final;
  void write(std::string, const std::string&) final;
  void write(std::string, const std::vector<double>&) final;
  void write(std::string, const std::vector<int>&) final;
  
  std::shared_ptr<sink_impl> group(std::string) final;
  
  ~h5file_sink_impl() = default;
};


class h5group_sink_impl : public sink_impl {
  const h5_group_write res;
  
  public:
  h5group_sink_impl(const h5_group_write res_);
  void write(std::string, const double&) final;
  void write(std::string, const int&) final;
  void write(std::string, const bool&) final;
  void write(std::string, const std::string&) final;
  void write(std::string, const std::vector<double>&) final;
  void write(std::string, const std::vector<int>&) final;
  
  std::shared_ptr<sink_impl> group(std::string) final;
  
  ~h5group_sink_impl() = default;
};




} //namespace detail



} // namespace EOS_Toolkit

#endif
