#ifndef INTERPOL_REGSPL_H
#define INTERPOL_REGSPL_H
#include <functional>
#include <memory>
#include <vector>
#include <array>

#include "config.h"
#include "interpol.h"
#include "datastore.h"

namespace EOS_Toolkit {
namespace detail {

class interpol_regspl_impl : public interpolator_impl {  
  public:

  struct segment {
    std::array<real_t, 4> c;
    segment(std::array<real_t, 4> c_) : c(c_) {};
    auto operator()(real_t t) const -> real_t;
    static auto hermite(real_t y0, real_t y1, real_t m0, real_t m1) 
    -> segment;
  };
  
  using interpolator_impl::func_t;
  using interpolator_impl::range_t;


  ///Default constructor. 
  interpol_regspl_impl()                    = default;
  ~interpol_regspl_impl()                   = default;
  interpol_regspl_impl(const interpol_regspl_impl&) = default;
  interpol_regspl_impl(interpol_regspl_impl&&);
  interpol_regspl_impl& operator=(interpol_regspl_impl);

  interpol_regspl_impl(std::vector<segment> segments, 
                       range_t range_x, range_t range_y);
  
  void swap(interpol_regspl_impl& other);
  
  ///Valid range. 
  auto range_x() const -> const range_t& final;

  ///Value range
  auto range_y() const -> const range_t& final;

  ///Look up value. 
  auto operator()(real_t x) const -> real_t final;
  
  void save(datasink s) const final;

  static auto from_vector(std::vector<real_t> values, range_t range_x) 
  -> interpol_regspl_impl;

  static auto from_function(func_t func, range_t range_x, 
                            size_t npoints)
  -> interpol_regspl_impl;

  static auto from_datasource(datasource src) 
  -> interpol_regspl_impl;

  auto transformed(func_t func) const
  -> interpol_regspl_impl;

  auto rescale_x(real_t scale) const
  -> interpol_regspl_impl;

  auto shift_x(real_t offset) const
  -> interpol_regspl_impl;

  auto make_rescale_x(real_t scale) const
  ->std::shared_ptr<interpolator_impl> final;
  
  auto make_transform(func_t) const
  ->std::shared_ptr<interpolator_impl> final;

  void assert_valid() const;

  static const std::string datastore_id;

  private:

  static auto get_dx(const range_t&, std::size_t) 
  -> real_t;
  
  static auto get_rgy(const std::vector<real_t>&)
  -> range_t;
  
  static auto make_seg(std::array<real_t, 4> y) 
  -> segment;

  std::vector<segment> segs;
  range_t rgx{0.,0.};
  range_t rgy{0.,0.};
  real_t dx{0};

};

void swap(interpol_regspl_impl& a, interpol_regspl_impl& b);

template<> struct source_proxy_reader<interpol_regspl_impl> {
  static void read(const datasource& s, std::string n, 
                   interpol_regspl_impl& t) 
  {
    t = interpol_regspl_impl::from_datasource(s / n);
  }
};


template<> struct source_proxy_writer<interpol_regspl_impl> {
  static void write(datasink& s, std::string n, 
                    const interpol_regspl_impl& t) 
  {
    t.save(s / n);
  }
};




} // namespace detail


auto make_interpol_regspl(detail::interpol_regspl_impl loc)
-> interpolator;

auto operator*(real_t, detail::interpol_regspl_impl)
->detail::interpol_regspl_impl;

auto operator*(detail::interpol_regspl_impl, real_t)
->detail::interpol_regspl_impl;

auto operator/(detail::interpol_regspl_impl, real_t)
->detail::interpol_regspl_impl;

auto operator/(real_t, detail::interpol_regspl_impl)
->detail::interpol_regspl_impl;



} // namespace EOS_Toolkit

#endif
