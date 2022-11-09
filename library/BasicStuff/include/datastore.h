#ifndef DATASTORE_H
#define DATASTORE_H

#include <memory>
#include <utility>
#include <string>
#include <vector>
#include <stdexcept>
#include <boost/optional.hpp>
#include "intervals.h"

namespace EOS_Toolkit {


namespace detail {

class source_impl {
  public:
  virtual bool has_data(std::string) const=0;
  virtual void read(std::string, double&) const=0;
  virtual void read(std::string, int&) const=0;
  virtual void read(std::string, bool&) const=0;
  virtual void read(std::string, std::string&) const=0;
  virtual void read(std::string, std::vector<double>&) const=0;
  virtual void read(std::string, std::vector<int>&) const=0;
  
  virtual bool has_group(std::string) const=0;
  virtual std::shared_ptr<source_impl> group(std::string) const=0;
  virtual ~source_impl() {}
};

class sink_impl {
  public:
  virtual void write(std::string, const double&) =0;
  virtual void write(std::string, const int&) =0;
  virtual void write(std::string, const bool&) =0;
  virtual void write(std::string, const std::string&) =0;
  virtual void write(std::string, const std::vector<double>&) =0;
  virtual void write(std::string, const std::vector<int>&) =0;
  
  virtual std::shared_ptr<sink_impl> group(std::string) =0;
  
  virtual ~sink_impl() {}
};


class source_proxy;
class sink_proxy;

} // namespace detail


class datasource {
  using impl_t = detail::source_impl;
  std::shared_ptr<impl_t> pimpl;
  std::shared_ptr<impl_t> parent;
  
  public:
  
  explicit datasource(std::shared_ptr<impl_t> pimpl_) 
  : pimpl{std::move(pimpl_)} 
  {}
  
  datasource(std::shared_ptr<impl_t> pimpl_, 
             std::shared_ptr<impl_t> parent_) 
  : pimpl{std::move(pimpl_)}, parent{std::move(parent_)}
  {}
  
  bool has_group(std::string n) const;
  datasource operator/(std::string n) const;
  bool has_data(std::string n) const;
  detail::source_proxy operator[](std::string n) const;
  
  template<class T> 
  void read(std::string n, T& t) const {
    pimpl->read(n, t);
  } 
};

class datasink {
  using impl_t = detail::sink_impl;
  std::shared_ptr<impl_t> pimpl;
  std::shared_ptr<impl_t> parent;
  
  public:
  
  explicit datasink(std::shared_ptr<impl_t> pimpl_) 
  : pimpl{pimpl_} 
  {}
   
  datasink(std::shared_ptr<impl_t> pimpl_, 
             std::shared_ptr<impl_t> parent_) 
  : pimpl{std::move(pimpl_)}, parent{std::move(parent_)}
  {}
  
  public:
  datasink operator/(std::string);
  detail::sink_proxy operator[](std::string);
  
  template<class T> 
  void write(std::string n, const T&t) {
    pimpl->write(n, t);
  } 
};

namespace detail {

template<class T> struct source_proxy_reader {
};

template<> struct source_proxy_reader<std::string> {
  static void read(const datasource& s, std::string n, std::string& t) {
    s.read(n, t);
  }
};
template<> struct source_proxy_reader<int> {
  static void read(const datasource& s, std::string n, int& t) {
    s.read(n, t);
  }
};
template<> struct source_proxy_reader<double> {
  static void read(const datasource& s, std::string n, double& t) {
    s.read(n, t);
  }
};
template<> struct source_proxy_reader<bool> {
  static void read(const datasource& s, std::string n, bool& t) {
    s.read(n, t);
  }
};
template<> struct source_proxy_reader<std::vector<double> > {
  static void read(const datasource& s, std::string n, 
                   std::vector<double>& t) {
    s.read(n, t);
  }
};
template<> struct source_proxy_reader<std::vector<int> > {
  static void read(const datasource& s, std::string n, 
                   std::vector<int>& t) {
    s.read(n, t);
  }
};

template<class T> struct source_proxy_reader<boost::optional<T> > {
  static void read(const datasource& s, std::string n, 
                   boost::optional<T>& o) {
    if (s.has_data(n)) {
      T t;
      source_proxy_reader<T>::read(s, n, t);
      o = std::move(t);
    }
    else {
      o = boost::none;
    }
  }
};




class source_proxy {
  const datasource& p;
  const std::string& n;
  
  public:
  source_proxy(const datasource& p_, const std::string& n_) 
  : p{p_}, n{n_} {}

  template<class T, 
          class R=decltype(&source_proxy_reader<T>::read) > 
  operator T() const {
    T t; 
    source_proxy_reader<T>::read(p, n, t); 
    return t;
  } 
};

template<class T> struct source_proxy_writer {
  static void write(datasink& s, std::string n, const T& t) {
    s.write(n, t);
  }
};


class sink_proxy {
  datasink& p;
  const std::string& n;
  
  public:
  sink_proxy(datasink& p_, const std::string& n_) 
  : p{p_}, n{n_} {}

  template<class T> 
  void operator=(const T& t) {
    source_proxy_writer<T>::write(p, n, t);
  }
  void operator=(const char* t) {
    source_proxy_writer<std::string>::write(p, n, std::string(t));
  }
};


template<> struct source_proxy_reader< interval<double> > {
  static void read(const datasource& s, std::string n, 
                   interval<double>& t) {
    auto g = s / n;
    t = interval<double>{ double(g["min"]), double(g["max"]) };
  }
};

template<> struct source_proxy_writer< interval<double> > {
  static void write(datasink& s, std::string n, 
                    const interval<double>& t) {
    auto g = s / n;
    g["min"] = t.min();
    g["max"] = t.max();
  }
};


} // namespace detail

} // namespace EOS_Toolkit


#endif



