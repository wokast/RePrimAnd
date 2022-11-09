#include "interpol_pchip_spline.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <iterator>
#include <stdexcept>

namespace {

template<class T>
bool is_strictly_increasing(const std::vector<T>& v)
{
  for (std::size_t i=1; i<v.size(); ++i) {
    if (v[i] <= v[i-1]) return false;
  }
  return true;
}

}

namespace EOS_Toolkit {
  
namespace detail {
  
const std::string interpol_pchip_impl::datastore_id {
  "pchip_spline"
};
  

interpol_pchip_impl::wrap_interp_accel::wrap_interp_accel()
: p{gsl_interp_accel_alloc()}
{
  if (p==nullptr) {
    throw std::runtime_error("interpol_pchip_impl: could not allocate memory");
  }
}

interpol_pchip_impl::wrap_interp_accel::~wrap_interp_accel()
{
  if (p!=nullptr) gsl_interp_accel_free(p);
}

interpol_pchip_impl::wrap_interp_cspline::wrap_interp_cspline( 
                   std::vector<double> x_, std::vector<double> y_)
: x{std::move(x_)}, y{std::move(y_)}
{
  const int min_points = 5;
  if (x.size() < min_points) {
    throw std::invalid_argument("interpol_pchip_impl: not enough "
                                "interpolation points");
  }
  if (x.size() != y.size()) {
    throw std::invalid_argument("interpol_pchip_impl: array size mismatch");
  }
  if (!is_strictly_increasing(x)) {
    throw std::runtime_error("interpol_pchip_impl: x-values must be strictly "
                             "increasing");
  }

  p = gsl_interp_alloc(gsl_interp_steffen, x.size());
  if (p == nullptr) {
    throw std::runtime_error("interpol_pchip_impl: could not allocate memory");
  }
  gsl_interp_init(p, &(x[0]), &(y[0]), x.size());
}

interpol_pchip_impl::wrap_interp_cspline::~wrap_interp_cspline()
{
  if (p!=nullptr) gsl_interp_free(p);
}
  
auto interpol_pchip_impl::wrap_interp_cspline::operator()(real_t t) 
const -> real_t
{
  assert(p);
  return gsl_interp_eval(p, &(x[0]), &(y[0]), t, acc.p);  
}


void interpol_pchip_impl::swap(interpol_pchip_impl& other)
{
  using std::swap;
  swap(spline, other.spline);
  swap(rgx, other.rgx);
  swap(rgy, other.rgy);
}

void swap(interpol_pchip_impl& a, interpol_pchip_impl& b)
{ 
  a.swap(b);
}


auto interpol_pchip_impl::get_rgx(const std::vector<real_t>& x)
-> range_t
{
  const int min_points = 5;
  if (x.size() < min_points) {
    throw std::invalid_argument("interpol_pchip_impl: not enough "
                                "sample points");
  }
  if (!is_strictly_increasing(x)) 
  {
    throw std::runtime_error("interpol_pchip_impl: sample positions " 
                             "must be strictly increasing");
  }
  return {x.front(), x.back()};
}

auto interpol_pchip_impl::get_rgy(const std::vector<real_t>& y)
-> range_t
{
  auto ext = std::minmax_element(y.begin(), y.end());
  return {*ext.first, *ext.second};
}



interpol_pchip_impl::interpol_pchip_impl(std::vector<real_t> sample_x, 
                                         std::vector<real_t> sample_y)
: rgx{get_rgx(sample_x)}, rgy{get_rgy(sample_y)},
  spline{
    std::make_shared<spline_t>(std::move(sample_x), 
                               std::move(sample_y))
  }
{}

void interpol_pchip_impl::assert_valid() const
{
  assert(spline);
}

auto interpol_pchip_impl::range_x() const -> const range_t& 
{
  assert_valid();
  return rgx;
}

auto interpol_pchip_impl::range_y() const -> const range_t&
{
  assert_valid();
  return rgy;
}

auto interpol_pchip_impl::from_vector(std::vector<real_t> sample_x, 
                                      std::vector<real_t> sample_y) 
-> interpol_pchip_impl
{
  return interpol_pchip_impl{std::move(sample_x), std::move(sample_y)};
}
  
auto interpol_pchip_impl::from_function(std::vector<real_t> x,
                                        func_t func)
-> interpol_pchip_impl
{
  std::vector<real_t> y;
  std::transform(x.begin(), x.end(), std::back_inserter(y),
                 func);
  
  return from_vector(std::move(x), std::move(y));
}

auto interpol_pchip_impl::transformed(func_t func) const
-> interpol_pchip_impl
{
  assert_valid();

  std::vector<real_t> yt;
  
  std::transform(spline->y.begin(), spline->y.end(), 
                 std::back_inserter(yt), func);
                 
  return from_vector(spline->x, std::move(yt));
}

auto interpol_pchip_impl::make_transform(func_t func) const
->std::shared_ptr<interpolator_impl>
{ 
  return std::make_shared<interpol_pchip_impl>(transformed(func));
}


auto interpol_pchip_impl::rescale_x(real_t scale) 
const -> interpol_pchip_impl
{       
  auto trans = [&] (real_t xold) {return xold*scale;};
  std::vector<real_t> xnew;
  std::transform(spline->x.begin(), spline->x.end(), 
                 std::back_inserter(xnew), trans);
  return interpol_pchip_impl{std::move(xnew), spline->y};
}

auto interpol_pchip_impl::make_rescale_x(real_t scale) const
->std::shared_ptr<interpolator_impl>
{
  return std::make_shared<interpol_pchip_impl>(rescale_x(scale));  
}


auto operator*(detail::interpol_pchip_impl i, real_t a)
->detail::interpol_pchip_impl
{
  return i.transformed([a](real_t y) {return a*y;});
}


auto operator*(real_t a, detail::interpol_pchip_impl i)
->detail::interpol_pchip_impl
{
  return i * a;
}

auto operator/(detail::interpol_pchip_impl i, real_t a)
->detail::interpol_pchip_impl
{
  return i * (1.0/a);
}

auto operator/(real_t a, detail::interpol_pchip_impl i)
->detail::interpol_pchip_impl
{
  return i.transformed([a](real_t y) {return a/y;});
}


auto interpol_pchip_impl::from_datasource(datasource s) 
-> interpol_pchip_impl
{
  std::string styp = s["interpolator_type"];
  if (styp != datastore_id) {
    throw std::runtime_error("unexpected interpolator type in "
                             "datasource encountered");
  }
  std::vector<real_t> x = s["sample_points"];
  std::vector<real_t> y = s["sample_values"];
  
  return from_vector(std::move(x), std::move(y));
}

void interpol_pchip_impl::save(datasink s) const
{
  assert_valid();
  s["interpolator_type"] = datastore_id;
  s["sample_points"] = spline->x;
  s["sample_values"] = spline->y;
}


/**
If x is outside the tabulated range, the function value at the 
closest boundary is returned
*/
real_t interpol_pchip_impl::operator()(real_t x) const
{
  assert_valid();
  
  return (*spline)(range_x().limit_to(x));
}



} // namespace detail


auto make_interpol_pchip_spline(
         const detail::interpol_pchip_impl& loc)
-> interpolator 
{
  return interpolator{
      std::make_shared<detail::interpol_pchip_impl>(loc)
  };  
}

auto make_interpol_pchip_spline(std::vector<real_t> x, 
                                std::vector<real_t> y)
-> interpolator 
{
  return make_interpol_pchip_spline(
          detail::interpol_pchip_impl::from_vector(std::move(x), 
                                                   std::move(y)));
}

auto make_interpol_pchip_spline(std::vector<real_t> x, 
                                std::function<real_t(real_t)> f)
-> interpolator 
{
  return make_interpol_pchip_spline(
          detail::interpol_pchip_impl::from_function(std::move(x), f));  
}

auto make_interpol_pchip_spline(datasource s)
-> interpolator 
{
  return make_interpol_pchip_spline(
          detail::interpol_pchip_impl::from_datasource(s));  
}



} // namespace EOS_Toolkit





