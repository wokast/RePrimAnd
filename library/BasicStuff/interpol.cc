#include "interpol_impl.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <stdexcept>

using namespace EOS_Toolkit;
using namespace EOS_Toolkit::detail;

namespace {

real_t get_log_map_offset(real_t a, real_t b, int mags)
{
  if (mags <= 0) {
    throw std::range_error("lookup_table_magx: magnitude bound "
                           "not strictly positive");
  }
  if (a<0) {
    throw std::range_error("lookup_table_magx: independent variable "
                           "range includes negative values");
  }

  real_t m10 = pow(10.0, -mags);
  real_t ofs = std::max(0.0, (m10*b - a) / (1.0 - m10));

  if (a+ofs <= 0) {
    throw std::range_error("lookup_table_magx: cannot handle "
                           "magnitude range");
  }

  return ofs;
}

template<class T>
bool is_strictly_increasing(const std::vector<T>& v)
{
  for (std::size_t i=1; i<v.size(); ++i) {
    if (v[i] <= v[i-1]) return false;
  }
  return true;
}

}

/**
Sample from arbitrary function over a given range with given
number of points.
**/
lookup_table::lookup_table(func_t func, range_t range, 
                           std::size_t npoints)
: y{}, rgx{range}
{
  if (npoints < 2) {
    throw std::range_error("lookup_table: need as least two "
                           "sample points");
  }
  real_t dx = range.length() / (npoints - 1.0);
  dxinv      = 1.0 / dx;

  for (std::size_t k=0; k<npoints; ++k) {
    real_t x = range.limit_to(range.min() + dx*k); 
    y.push_back(func(x));
  }
  
  auto ext = std::minmax_element(y.begin(), y.end());
  rgy      = {*ext.first, *ext.second};

}

/**
If x is outside the tabulated range, the function value at the 
closest boundary is returned
*/
real_t lookup_table::operator()(real_t x) const
{
  x = range_x().limit_to(x);
  const real_t s  = (x - range_x().min()) * dxinv;
  
  assert(s >= 0);
  const unsigned int i = floor(s);
  
  const unsigned int j = i + 1;
  if (j >= y.size()) {  //can happen only by rounding errors
    return y.back();
  }
  return  (s-i) * y[j] + (j-s) * y[i];
}



lookup_table_magx::lookup_table_magx(func_t func, range_t range, 
                                     size_t npoints, int magnitudes)
: rgx{range},
  x_offs{get_log_map_offset(range.min(), range.max(), magnitudes)}
{
  auto gunc = [this, &func] (real_t lgx) {
    return func(exp(lgx) - x_offs);
  };

  range_t lgrg{log(rgx.min() + x_offs), 
               log(rgx.max() + x_offs)};
  
  tbl = {gunc, lgrg, npoints};
}

/**
If x is outside the tabulated range, the function value at the 
closest boundary is returned
*/
real_t lookup_table_magx::operator()(real_t x) const
{
  x = range_x().limit_to(x);
  return tbl(log(x + x_offs));
}



cspline_mono::cspline_mono(const std::vector<real_t>& x_, 
                           const std::vector<real_t>& y_)
{
  std::vector<double> x;
  std::vector<double> y;
  std::copy(x_.begin(), x_.end(), std::back_inserter(x));
  std::copy(y_.begin(), y_.end(), std::back_inserter(y));
    
  pimpl = std::make_shared<cspline_mono_impl>(
                         std::move(x), std::move(y));  
  
}

auto cspline_mono::valid() const 
-> const cspline_mono_impl&
{
  if (!pimpl) {
    throw std::logic_error("cspline_mono: uninitialized use.");
  }
  return *pimpl;
}

auto cspline_mono::range_x() const -> const range_t& 
{
  return valid().range_x();
}

auto cspline_mono::range_y() const -> const range_t& 
{
  return valid().range_y();
}

real_t cspline_mono::operator()(real_t x) const
{
  return valid()(x);
}


wrap_interp_accel::wrap_interp_accel()
: p{gsl_interp_accel_alloc()}
{
  if (p==nullptr) {
    throw std::runtime_error("cspline_mono: could not allocate memory");
  }
}

wrap_interp_accel::~wrap_interp_accel()
{
  gsl_interp_accel_free(p);
}

detail::wrap_interp_cspline::wrap_interp_cspline( 
  const std::vector<double>& x, const std::vector<double>& y)
{
  const int min_points = 5;
  if (x.size() < min_points) {
    throw std::invalid_argument("cspline_mono: not enough "
                                "interpolation points");
  }
  if (x.size() != y.size()) {
    throw std::invalid_argument("cspline_mono: array size mismatch");
  }
  if (!is_strictly_increasing(x)) {
    throw std::runtime_error("cspline_mono: x-values must be strictly "
                             "increasing");
  }

  p = gsl_interp_alloc(gsl_interp_steffen, x.size());
  if (p == nullptr) {
    throw std::runtime_error("cspline_mono: could not allocate memory");
  }
  gsl_interp_init(p, &(x[0]), &(y[0]), x.size());
}

wrap_interp_cspline::~wrap_interp_cspline()
{
  gsl_interp_free(p);
}

detail::cspline_mono_impl::cspline_mono_impl(
                std::vector<double> x_, std::vector<double> y_)
: x{std::move(x_)}, y{std::move(y_)}, interp{x,y}
{  
  auto ext = std::minmax_element(y.begin(), y.end());
  rgy = {*ext.first, *ext.second};
  rgx = {x.front(), x.back()};
}

real_t cspline_mono_impl::operator()(real_t t) const
{
  t = range_x().limit_to(t);
  return gsl_interp_eval(interp.p, &(x[0]), &(y[0]), t, acc.p);
}









