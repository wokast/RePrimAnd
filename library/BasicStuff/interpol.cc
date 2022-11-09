#include "interpol.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <stdexcept>

namespace EOS_Toolkit {
//~ namespace detail {
  
interpolator::interpolator(std::shared_ptr<interpolator_impl> pimpl_)
: pimpl{pimpl_}
{}

auto interpolator::valid() const 
-> const interpolator_impl&
{
  if (!pimpl) {
    throw std::logic_error("interpolator: uninitialized use.");
  }
  return *pimpl;
}

auto interpolator::range_x() const -> const range_t& 
{
  return valid().range_x();
}

auto interpolator::range_y() const ->const range_t&
{
  return valid().range_y();
}


auto interpolator::transformed(func_t f) const ->interpolator
{
  return interpolator{make_transform(f)};
}

auto interpolator::make_transform(func_t f) const
->std::shared_ptr<interpolator_impl>
{
  return valid().make_transform(f);
}

auto interpolator::rescale_x(real_t scale) const
-> interpolator
{
  return interpolator{make_rescale_x(scale)};
}

auto interpolator::make_rescale_x(real_t scale) const
->std::shared_ptr<interpolator_impl>
{
  return valid().make_rescale_x(scale);  
}

  
//~ }


auto operator*(real_t a, interpolator i)
->interpolator
{
  return i.transformed([a](real_t y) {return a*y;});
}

auto operator*(interpolator i, real_t a)
->interpolator
{
  return a*i;
}

auto operator/(interpolator i, real_t a)
->interpolator
{
  return i * (1.0/a);
}

auto operator/(real_t a, interpolator i)
->interpolator
{
  return i.transformed([a](real_t y) {return a/y;});
}

}



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





