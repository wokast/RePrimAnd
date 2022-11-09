#ifndef INTERPOL_LINEAR_H
#define INTERPOL_LINEAR_H
#include <functional>
#include <memory>
#include <vector>
#include "config.h"
#include "interpol.h"
#include "datastore.h"

namespace EOS_Toolkit {


namespace detail {

class interpol_reglin_impl : public interpolator_impl {  
  public:

  using interpolator_impl::func_t;
  using interpolator_impl::range_t;


  ///Default constructor. 
  interpol_reglin_impl()                    = default;
  ~interpol_reglin_impl()                   = default;
  interpol_reglin_impl(const interpol_reglin_impl&) = default;
  interpol_reglin_impl(interpol_reglin_impl&&);
  interpol_reglin_impl& operator=(interpol_reglin_impl);

  interpol_reglin_impl(std::vector<real_t> values, range_t range_x);
  
  void swap(interpol_reglin_impl& other);
  
  ///Valid range. 
  auto range_x() const -> const range_t& final;

  ///Value range
  auto range_y() const -> const range_t& final;

  ///Look up value. 
  auto operator()(real_t x) const -> real_t final;
  
  void save(datasink s) const final;

  static auto from_vector(std::vector<real_t> values, range_t range_x) 
  -> interpol_reglin_impl;

  static auto from_function(func_t func, range_t range_x, 
                            size_t npoints)
  -> interpol_reglin_impl;

  static auto from_datasource(datasource src) 
  -> interpol_reglin_impl;

  auto transformed(func_t func) const
  -> interpol_reglin_impl;
  
  auto make_transform(func_t) const
  ->std::shared_ptr<interpolator_impl> final;

  auto rescale_x(real_t scale) const
  -> interpol_reglin_impl;

  auto shift_x(real_t offset) const
  -> interpol_reglin_impl;

  auto make_rescale_x(real_t scale) const
  ->std::shared_ptr<interpolator_impl> final;

  void assert_valid() const;
 
  static const std::string datastore_id;
 
  private:

  static auto get_dx(const range_t&, std::size_t) 
  -> real_t;
  
  static auto get_rgy(const std::vector<real_t>&)
  -> range_t;

  
  std::vector<real_t> y;
  real_t dxinv{0.0};
  range_t rgx{0,0};
  range_t rgy{0,0};
};

void swap(interpol_reglin_impl& a, interpol_reglin_impl& b);

template<> struct source_proxy_reader<interpol_reglin_impl> {
  static void read(const datasource& s, std::string n, 
                   interpol_reglin_impl& t) 
  {
    t = interpol_reglin_impl::from_datasource(s / n);
  }
};


template<> struct source_proxy_writer<interpol_reglin_impl> {
  static void write(datasink& s, std::string n, 
                    const interpol_reglin_impl& t) 
  {
    t.save(s / n);
  }
};


class interpol_loglin_impl : public interpolator_impl {
  public:

  using func_t  = std::function<real_t(real_t)>;
  using interpolator_impl::range_t;


  ///Default constructor. 
  interpol_loglin_impl()                    = default;
  ~interpol_loglin_impl()                   = default;
  
  interpol_loglin_impl(const interpol_loglin_impl&) = default;
  interpol_loglin_impl(interpol_loglin_impl&&) = default;
  
  auto operator=(const interpol_loglin_impl&) 
  -> interpol_loglin_impl& = default;
  
  auto operator=(interpol_loglin_impl&&) 
  ->interpol_loglin_impl& = default;

  interpol_loglin_impl(interpol_reglin_impl yz_);

  void swap(interpol_loglin_impl& other);
  
  ///Valid range. 
  auto range_x() const -> const range_t& final;

  ///Value range
  auto range_y() const -> const range_t& final;

  ///Look up value. 
  auto operator()(real_t x) const -> real_t final;
  
  static auto from_vector(std::vector<real_t> values, range_t range_x) 
  -> interpol_loglin_impl;

  static auto from_function(func_t func, range_t range_x, 
                            size_t npoints)
  -> interpol_loglin_impl;
  
  void save(datasink dst) const final;

  static auto from_datasource(datasource src) 
  -> interpol_loglin_impl;

  auto transformed(func_t func) const
  -> interpol_loglin_impl;
  
  auto make_transform(func_t) const
  ->std::shared_ptr<interpolator_impl> final;

  auto rescale_x(real_t scale) const
  -> interpol_loglin_impl;

  auto make_rescale_x(real_t scale) const
  ->std::shared_ptr<interpolator_impl> final;

  static const std::string datastore_id;

  private:

  void assert_valid() const;

  static auto x2z(real_t x) -> real_t;
  static auto z2x(real_t) -> real_t;
  
  static auto rgz2rgx(range_t rgz) -> range_t;
  static auto rgx2rgz(range_t rgx) -> range_t;
  
  //~ static auto get_log_map_offset(real_t a, real_t b, int mags) 
  //~ ->real_t;
  
  
  //~ real_t x_offs{1.0};
  interpol_reglin_impl yz;  
  range_t rgx{0,0};
  
};

void swap(interpol_loglin_impl& a, interpol_loglin_impl& b);

template<> struct source_proxy_reader<interpol_loglin_impl> {
  static void read(const datasource& s, std::string n, 
                   interpol_loglin_impl& t) 
  {
    t = interpol_loglin_impl::from_datasource(s / n);
  }
};


template<> struct source_proxy_writer<interpol_loglin_impl> {
  static void write(datasink& s, std::string n, 
                    const interpol_loglin_impl& t) 
  {
    t.save(s / n);
  }
};

} // namespace detail


auto operator*(real_t, detail::interpol_reglin_impl)
->detail::interpol_reglin_impl;

auto operator*(detail::interpol_reglin_impl, real_t)
->detail::interpol_reglin_impl;

auto operator/(detail::interpol_reglin_impl, real_t)
->detail::interpol_reglin_impl;

auto operator/(real_t, detail::interpol_reglin_impl)
->detail::interpol_reglin_impl;

auto operator*(real_t, detail::interpol_loglin_impl)
->detail::interpol_loglin_impl;

auto operator*(detail::interpol_loglin_impl, real_t)
->detail::interpol_loglin_impl;

auto operator/(detail::interpol_loglin_impl, real_t)
->detail::interpol_loglin_impl;

auto operator/(real_t, detail::interpol_loglin_impl)
->detail::interpol_loglin_impl;

} // namespace EOS_Toolkit

#endif
