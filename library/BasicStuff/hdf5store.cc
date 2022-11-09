#include "hdf5store.h"
#include "hdf5store_impl.h"
#include <cassert>
#include <stdexcept>


namespace EOS_Toolkit {

namespace detail {



h5file_source_impl::h5file_source_impl(std::string path)
: res(path) 
{}

bool h5file_source_impl::has_data(std::string n) const
{
  return has_entry(res, n);
}

void h5file_source_impl::read(std::string n, double& v) const
{
  read_attr(res, n, v);
}

void h5file_source_impl::read(std::string n, int& v) const
{
  read_attr(res, n, v);
}

void h5file_source_impl::read(std::string n, bool& v) const
{
  read_attr(res, n, v);
}

void h5file_source_impl::read(std::string n, std::string& v) const
{
  read_attr(res, n, v);
}

void h5file_source_impl::read(std::string n, 
                              std::vector<double>& v) const 
{
  read_data(res, n, v);
}

void h5file_source_impl::read(std::string n, 
                              std::vector<int>& v) const
{
  read_data(res, n, v);
}


bool h5file_source_impl::has_group(std::string n) const
{
  return has_entry(res, n); // Todo: check specifically for groups 
}

std::shared_ptr<source_impl> 
h5file_source_impl::group(std::string n) const 
{
  return std::make_shared<detail::h5group_source_impl>(
                                        h5_group_read(res, n));
}



h5group_source_impl::h5group_source_impl(const h5_group_read res_)
: res{res_}
{}


bool h5group_source_impl::has_data(std::string n) const
{
  return has_entry(res, n);
}

void h5group_source_impl::read(std::string n, double& v) const
{
  read_attr(res, n, v);
}

void h5group_source_impl::read(std::string n, int& v) const
{
  read_attr(res, n, v);
}

void h5group_source_impl::read(std::string n, bool& v) const
{
  read_attr(res, n, v);
}

void h5group_source_impl::read(std::string n, std::string& v) const
{
  read_attr(res, n, v);
}

void h5group_source_impl::read(std::string n, 
                              std::vector<double>& v) const 
{
  read_data(res, n, v);
}

void h5group_source_impl::read(std::string n, 
                              std::vector<int>& v) const
{
  read_data(res, n, v);
}


bool h5group_source_impl::has_group(std::string n) const
{
  return has_entry(res, n); // Todo: check specifically for groups 
}


std::shared_ptr<source_impl> 
h5group_source_impl::group(std::string n) const 
{
  return std::make_shared<detail::h5group_source_impl>(
                                          h5_group_read(res, n));
}








h5file_sink_impl::h5file_sink_impl(std::string path)
: res(path) 
{}


void h5file_sink_impl::write(std::string n, const double& v) 
{
  write_attr(res, n, v);
}

void h5file_sink_impl::write(std::string n, const int& v) 
{
  write_attr(res, n, v);
}

void h5file_sink_impl::write(std::string n, const bool& v) 
{
  write_attr(res, n, v);
}

void h5file_sink_impl::write(std::string n, const std::string& v) 
{
  write_attr(res, n, v);
}

void h5file_sink_impl::write(std::string n, 
                             const std::vector<double>& v)
{
  write_data(res, n, v);
}

void h5file_sink_impl::write(std::string n, 
                             const std::vector<int>& v) 
{
  write_data(res, n, v);
}



std::shared_ptr<sink_impl> 
h5file_sink_impl::group(std::string n) 
{
  return std::make_shared<detail::h5group_sink_impl>(
                                          h5_group_write(res, n));
}





h5group_sink_impl::h5group_sink_impl(const h5_group_write res_)
: res(res_) 
{}


void h5group_sink_impl::write(std::string n, const double& v) 
{
  write_attr(res, n, v);
}

void h5group_sink_impl::write(std::string n, const int& v) 
{
  write_attr(res, n, v);
}

void h5group_sink_impl::write(std::string n, const bool& v) 
{
  write_attr(res, n, v);
}

void h5group_sink_impl::write(std::string n, const std::string& v) 
{
  write_attr(res, n, v);
}

void h5group_sink_impl::write(std::string n, 
                              const std::vector<double>& v)
{
  write_data(res, n, v);
}

void h5group_sink_impl::write(std::string n, 
                              const std::vector<int>& v) 
{
  write_data(res, n, v);
}


std::shared_ptr<sink_impl> 
h5group_sink_impl::group(std::string n) 
{
  return std::make_shared<detail::h5group_sink_impl>(
                                          h5_group_write(res, n));
}






} // namespace detail 


datasource make_hdf5_file_source(std::string path)
{
  return datasource(std::make_shared<detail::h5file_source_impl>(path));
}

datasink make_hdf5_file_sink(std::string path)
{
  return datasink(std::make_shared<detail::h5file_sink_impl>(path));
}

} // namespace EOS_Toolkit 

