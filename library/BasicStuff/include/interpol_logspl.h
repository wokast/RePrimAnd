#ifndef INTERPOL_LOGSPL_H
#define INTERPOL_LOGSPL_H
#include <functional>
#include <memory>
#include <vector>

#include "config.h"
#include "interpol.h"
#include "interpol_regspl.h"
#include "datastore.h"

namespace EOS_Toolkit {
namespace detail {


class interpol_logspl_impl : public interpolator_impl {
  public:

  using func_t  = std::function<real_t(real_t)>;
  using interpolator_impl::range_t;


  ///Default constructor. 
  interpol_logspl_impl()                    = default;
  ~interpol_logspl_impl()                   = default;
  
  interpol_logspl_impl(const interpol_logspl_impl&) = default;
  interpol_logspl_impl(interpol_logspl_impl&&) = default;
  
  auto operator=(const interpol_logspl_impl&) 
  -> interpol_logspl_impl& = default;
  
  auto operator=(interpol_logspl_impl&&) 
  ->interpol_logspl_impl& = default;

  interpol_logspl_impl(interpol_regspl_impl yz_);

  void swap(interpol_logspl_impl& other);
  
  ///Valid range. 
  auto range_x() const -> const range_t& final;

  ///Value range
  auto range_y() const -> const range_t& final;

  ///Look up value. 
  auto operator()(real_t x) const -> real_t final;
  
  static auto from_vector(std::vector<real_t> values, range_t range_x) 
  -> interpol_logspl_impl;

  static auto from_function(func_t func, range_t range_x, 
                            size_t npoints)
  -> interpol_logspl_impl;
  
  void save(datasink dst) const final;

  static auto from_datasource(datasource src) 
  -> interpol_logspl_impl;

  auto transformed(func_t func) const
  -> interpol_logspl_impl;

  auto make_transform(func_t) const
  ->std::shared_ptr<interpolator_impl> final;

  auto rescale_x(real_t scale) const
  -> interpol_logspl_impl;

  auto shift_x(real_t offset) const
  -> interpol_logspl_impl;
  
  auto make_rescale_x(real_t scale) const
  ->std::shared_ptr<interpolator_impl> final;
  

  static const std::string datastore_id;

  static auto x2z(real_t) -> real_t;
  static auto z2x(real_t) -> real_t;
  
  static auto rgz2rgx(range_t) -> range_t;
  static auto rgx2rgz(range_t) -> range_t;
  
  //~ static auto get_log_map_offset(real_t a, real_t b, int mags) 
  //~ ->real_t;

  void assert_valid() const;

  private:
  
  interpol_regspl_impl yz;  
  range_t rgx{0,0};
  
};

void swap(interpol_logspl_impl& a, interpol_logspl_impl& b);

template<> struct source_proxy_reader<interpol_logspl_impl> {
  static void read(const datasource& s, std::string n, 
                   interpol_logspl_impl& t) 
  {
    t = interpol_logspl_impl::from_datasource(s / n);
  }
};


template<> struct source_proxy_writer<interpol_logspl_impl> {
  static void write(datasink& s, std::string n, 
                    const interpol_logspl_impl& t) 
  {
    t.save(s / n);
  }
};



class interpol_llogspl_impl : public interpolator_impl {
  public:

  using func_t  = std::function<real_t(real_t)>;
  using interpolator_impl::range_t;


  ///Default constructor. 
  interpol_llogspl_impl()                    = default;
  ~interpol_llogspl_impl()                   = default;
  
  interpol_llogspl_impl(const interpol_llogspl_impl&) = default;
  interpol_llogspl_impl(interpol_llogspl_impl&&) = default;
  
  auto operator=(const interpol_llogspl_impl&) 
  -> interpol_llogspl_impl& = default;
  
  auto operator=(interpol_llogspl_impl&&) 
  ->interpol_llogspl_impl& = default;

  interpol_llogspl_impl(interpol_logspl_impl yz_);

  void swap(interpol_llogspl_impl& other);
  
  ///Valid range. 
  auto range_x() const -> const range_t& final;

  ///Value range
  auto range_y() const -> const range_t& final;

  ///Look up value. 
  auto operator()(real_t x) const -> real_t final;
  
  static auto from_vector(std::vector<real_t> values, range_t range_x) 
  -> interpol_llogspl_impl;

  static auto from_function(func_t func, range_t range_x, 
                            size_t npoints)
  -> interpol_llogspl_impl;
  
  void save(datasink dst) const final;

  static auto from_datasource(datasource src) 
  -> interpol_llogspl_impl;

  auto transformed(func_t func) const
  -> interpol_llogspl_impl;

  auto make_transform(func_t) const
  ->std::shared_ptr<interpolator_impl> final;

  auto rescale_x(real_t scale) const
  -> interpol_llogspl_impl;
  
  auto make_rescale_x(real_t scale) const
  ->std::shared_ptr<interpolator_impl> final;
  
  static const std::string datastore_id;

  private:

  void assert_valid() const;
  
  interpol_logspl_impl yz;  
  range_t rgy{0,0};
  
};

void swap(interpol_llogspl_impl& a, interpol_llogspl_impl& b);

template<> struct source_proxy_reader<interpol_llogspl_impl> {
  static void read(const datasource& s, std::string n, 
                   interpol_llogspl_impl& t) 
  {
    t = interpol_llogspl_impl::from_datasource(s / n);
  }
};


template<> struct source_proxy_writer<interpol_llogspl_impl> {
  static void write(datasink& s, std::string n, 
                    const interpol_llogspl_impl& t) 
  {
    t.save(s / n);
  }
};

} // namespace detail



auto make_interpol_logspl(detail::interpol_logspl_impl loc)
-> interpolator; 

auto operator*(real_t, detail::interpol_logspl_impl)
->detail::interpol_logspl_impl;

auto operator*(detail::interpol_logspl_impl, real_t)
->detail::interpol_logspl_impl;

auto operator/(detail::interpol_logspl_impl, real_t)
->detail::interpol_logspl_impl;

auto operator/(real_t, detail::interpol_logspl_impl)
->detail::interpol_logspl_impl;




auto make_interpol_llogspl(detail::interpol_llogspl_impl loc)
-> interpolator; 

auto operator*(real_t, detail::interpol_llogspl_impl)
->detail::interpol_llogspl_impl;

auto operator*(detail::interpol_llogspl_impl, real_t)
->detail::interpol_llogspl_impl;

auto operator/(detail::interpol_llogspl_impl, real_t)
->detail::interpol_llogspl_impl;

auto operator/(real_t, detail::interpol_llogspl_impl)
->detail::interpol_llogspl_impl;

} // namespace EOS_Toolkit

#endif
