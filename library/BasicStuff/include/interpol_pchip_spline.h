#ifndef INTERPOL_PCHIP_SPLINE_H
#define INTERPOL_PCHIP_SPLINE_H
#include <functional>
#include <memory>
#include <vector>
#include <gsl/gsl_spline.h>
#include "config.h"
#include "interpol.h"
#include "datastore.h"

namespace EOS_Toolkit {


namespace detail {




class interpol_pchip_impl : public interpolator_impl {  
      
  struct wrap_interp_accel {
    gsl_interp_accel* p{nullptr};
    
    wrap_interp_accel();
    ~wrap_interp_accel();
  };

  struct wrap_interp_cspline {
    wrap_interp_accel acc;
    gsl_interp* p{nullptr};
    std::vector<real_t> x;
    std::vector<real_t> y;
    
    wrap_interp_cspline() = delete;
    wrap_interp_cspline(const wrap_interp_cspline&) = delete;
    wrap_interp_cspline(wrap_interp_cspline&&) = delete;
    wrap_interp_cspline& operator=(const wrap_interp_cspline&) = delete;
    wrap_interp_cspline& operator=(wrap_interp_cspline&&) = delete;
    wrap_interp_cspline(std::vector<double> x_, 
                        std::vector<double> y_);
    auto operator()(real_t t) const -> real_t;
    
    ~wrap_interp_cspline();
  };
    
  public:
  
  using spline_t = wrap_interp_cspline; 
  using interpolator_impl::func_t;
  using interpolator_impl::range_t;

  ///Default constructor. 
  interpol_pchip_impl()                    = default;
  ~interpol_pchip_impl()                   = default;
  interpol_pchip_impl(const interpol_pchip_impl&) = default;
  interpol_pchip_impl(interpol_pchip_impl&&) = default;
  interpol_pchip_impl& operator=(const interpol_pchip_impl&) = default;
  interpol_pchip_impl& operator=(interpol_pchip_impl&&) = default;

  interpol_pchip_impl(std::vector<real_t> sample_x, 
                      std::vector<real_t> sample_y);
  
  void swap(interpol_pchip_impl& other);
  
  ///Valid range. 
  auto range_x() const -> const range_t& final;

  ///Value range
  auto range_y() const -> const range_t& final;

  ///Look up value. 
  auto operator()(real_t x) const -> real_t final;
  
  void save(datasink s) const;

  static auto get_rgx(const std::vector<real_t>&) -> range_t;
  static auto get_rgy(const std::vector<real_t>&) -> range_t;
  
  static auto from_vector(std::vector<real_t> sample_x, 
                          std::vector<real_t> sample_y) 
  -> interpol_pchip_impl;

  static auto from_function(std::vector<real_t> sample_x, func_t func)
  -> interpol_pchip_impl;

  static auto from_datasource(datasource src) 
  -> interpol_pchip_impl;

  auto transformed(func_t func) const
  -> interpol_pchip_impl;
  
  auto make_transform(func_t) const
  -> std::shared_ptr<interpolator_impl> final;

  auto rescale_x(real_t scale) const
  -> interpol_pchip_impl;

  auto make_rescale_x(real_t scale) const
  -> std::shared_ptr<interpolator_impl> final;
  
  static const std::string datastore_id;

  private:

  void assert_valid() const;

  range_t rgx;
  range_t rgy;

  std::shared_ptr<spline_t> spline;
};

void swap(interpol_pchip_impl& a, interpol_pchip_impl& b);

template<> struct source_proxy_reader<interpol_pchip_impl> {
  static void read(const datasource& s, std::string n, 
                   interpol_pchip_impl& t) 
  {
    t = interpol_pchip_impl::from_datasource(s / n);
  }
};


template<> struct source_proxy_writer<interpol_pchip_impl> {
  static void write(datasink& s, std::string n, 
                    const interpol_pchip_impl& t) 
  {
    t.save(s / n);
  }
};


} // namespace detail


auto operator*(real_t, detail::interpol_pchip_impl)
->detail::interpol_pchip_impl;

auto operator*(detail::interpol_pchip_impl, real_t)
->detail::interpol_pchip_impl;

auto operator/(detail::interpol_pchip_impl, real_t)
->detail::interpol_pchip_impl;

auto operator/(real_t, detail::interpol_pchip_impl)
->detail::interpol_pchip_impl;

} // namespace EOS_Toolkit

#endif
