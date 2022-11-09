#include "interpol_linear.h"
#include <cmath>
#include <algorithm>
#include <iterator>
#include <stdexcept>



namespace EOS_Toolkit {
  
namespace detail {

const std::string interpol_reglin_impl::datastore_id {
  "regular_spaced_linear"
};

interpol_reglin_impl::interpol_reglin_impl(
                             interpol_reglin_impl&& other)
{
  other.swap(*this);
}

interpol_reglin_impl& interpol_reglin_impl::operator=(
                                        interpol_reglin_impl other)
{
  other.swap(*this);
  return *this;
}

void interpol_reglin_impl::swap(interpol_reglin_impl& other)
{
  using std::swap;
  swap(y, other.y);
  swap(dxinv, other.dxinv);
  swap(rgx, other.rgx);
  swap(rgy, other.rgy);
}

void swap(interpol_reglin_impl& a, interpol_reglin_impl& b)
{ 
  a.swap(b);
}

void interpol_reglin_impl::assert_valid() const
{
  assert(!y.empty());
}

auto interpol_reglin_impl::get_dx(const range_t& rgx, 
                                     std::size_t npoints)
-> real_t
{
  if (npoints < 2) {
    throw std::range_error("interpol_reglin_impl: need as least two "
                           "sample points");
  }
  
  real_t dx = rgx.length() / (npoints - 1.0);
  
  if (dx <= 0) {
    throw std::range_error("interpol_reglin_impl: degenerate x-range");
  }
  
  return dx;
}

auto interpol_reglin_impl::get_rgy(const std::vector<real_t>& y)
-> range_t
{
  auto ext = std::minmax_element(y.begin(), y.end());
  return {*ext.first, *ext.second};
}

interpol_reglin_impl::interpol_reglin_impl(std::vector<real_t> values,
                                           range_t range_x)
: y{values}, dxinv{1.0/get_dx(range_x, values.size())}, 
  rgx{range_x}, rgy{get_rgy(values)}
{}

auto interpol_reglin_impl::range_x() const -> const range_t& 
{
  assert_valid();
  return rgx;
}

auto interpol_reglin_impl::range_y() const -> const range_t&
{
  assert_valid();
  return rgy;
}

auto interpol_reglin_impl::from_vector(std::vector<real_t> values, 
                                       range_t range_x) 
-> interpol_reglin_impl
{
  return interpol_reglin_impl{std::move(values), range_x};
}
  
auto interpol_reglin_impl::from_function(func_t func, range_t range_x, 
                                         size_t npoints)
-> interpol_reglin_impl
{
  real_t dx = get_dx(range_x, npoints);
  std::vector<real_t> y;
  
  for (std::size_t k=0; k<npoints; ++k) {
    real_t x = range_x.limit_to(range_x.min() + dx*k); 
    y.push_back(func(x));
  }
  return from_vector(std::move(y), range_x);
}

auto interpol_reglin_impl::transformed(func_t func) const
-> interpol_reglin_impl
{
  assert_valid();

  std::vector<real_t> yt;
  
  std::transform(y.begin(), y.end(), std::back_inserter(yt),
                 func);
                 
  return from_vector(std::move(yt), rgx);
}

auto interpol_reglin_impl::make_transform(func_t func) const
->std::shared_ptr<interpolator_impl>
{ 
  return std::make_shared<interpol_reglin_impl>(transformed(func));
}

auto interpol_reglin_impl::rescale_x(real_t scale) 
const -> interpol_reglin_impl
{
  auto xnew = [&] (real_t xold) {return xold*scale;};
  range_t rgxnew{ xnew(rgx.min()), xnew(rgx.max()) };
  return from_vector(y, rgxnew);
}

auto interpol_reglin_impl::shift_x(real_t offset) 
const -> interpol_reglin_impl
{
  auto xnew = [&] (real_t xold) {return xold + offset;};
  range_t rgxnew{ xnew(rgx.min()), xnew(rgx.max()) };
  return from_vector(y, rgxnew);
}

auto interpol_reglin_impl::make_rescale_x(real_t scale) const
->std::shared_ptr<interpolator_impl>
{
  return std::make_shared<interpol_reglin_impl>(rescale_x(scale));  
}


auto interpol_reglin_impl::from_datasource(datasource s) 
-> interpol_reglin_impl
{
  std::string styp = s["interpolator_type"];
  if (styp != datastore_id) {
    throw std::runtime_error("unexpected interpolator type in "
                             "datasource encountered");
  }
  std::vector<real_t> y = s["sample_values"];
  interval<real_t> rg{(real_t)s["range_min"], (real_t)s["range_max"]};
  
  return {std::move(y), rg};
}

void interpol_reglin_impl::save(datasink s) const
{
  assert_valid();
  
  s["interpolator_type"] = datastore_id;
  s["sample_values"] = y;
  s["range_min"] = rgx.min();
  s["range_max"] = rgx.max();
}


/**
If x is outside the tabulated range, the function value at the 
closest boundary is returned
*/
real_t interpol_reglin_impl::operator()(real_t x) const
{
  assert_valid();
  
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



const std::string interpol_loglin_impl::datastore_id {
  "log_spaced_linear"
};

auto interpol_loglin_impl::x2z(real_t x)
-> real_t
{
  return log(x);
}

auto interpol_loglin_impl::z2x(real_t z)
-> real_t
{
  return exp(z);
}

auto interpol_loglin_impl::rgz2rgx(range_t rgz)
-> range_t
{
  return {z2x(rgz.min()), z2x(rgz.max())};
  
}

auto interpol_loglin_impl::rgx2rgz(range_t rgx)
-> range_t
{
  return {x2z(rgx.min()), x2z(rgx.max())};
}


//~ auto interpol_loglin_impl::get_log_map_offset(real_t a, 
                                              //~ real_t b, int mags)
//~ -> real_t 
//~ {
  //~ if (mags <= 0) {
    //~ throw std::range_error("interpol_loglin_impl: magnitude bound "
                           //~ "not strictly positive");
  //~ }
  //~ if (a<0) {
    //~ throw std::range_error("interpol_loglin_impl: independent variable "
                           //~ "range includes negative values");
  //~ }

  //~ real_t m10 = pow(10.0, -mags);
  //~ real_t ofs = std::max(0.0, (m10*b - a) / (1.0 - m10));

  //~ if (a+ofs <= 0) {
    //~ throw std::range_error("interpol_loglin_impl: cannot handle "
                           //~ "magnitude range");
  //~ }

  //~ return ofs;
//~ }



interpol_loglin_impl::interpol_loglin_impl(interpol_reglin_impl yz_)
: yz{std::move(yz_)}, rgx{rgz2rgx(yz.range_x())}
{}


void interpol_loglin_impl::swap(interpol_loglin_impl& other)
{
  using std::swap;
  swap(yz, other.yz);
  swap(rgx, other.rgx);
}

void swap(interpol_loglin_impl& a, interpol_loglin_impl& b)
{ 
  a.swap(b);
}


void interpol_loglin_impl::assert_valid() const
{
  yz.assert_valid();
}

auto interpol_loglin_impl::range_x() const -> const range_t& 
{
  assert_valid();
  return rgx;
}

auto interpol_loglin_impl::range_y() const -> const range_t&
{
  return yz.range_y();
}


auto interpol_loglin_impl::operator()(real_t x) const -> real_t
{
  return yz(x2z(x));
}

auto interpol_loglin_impl::from_vector(std::vector<real_t> values, 
                                       range_t range_x) 
  -> interpol_loglin_impl
{
  auto rgz = rgx2rgz(range_x);
  interpol_reglin_impl yz{std::move(values), rgz};
  return interpol_loglin_impl{std::move(yz)};
}

auto interpol_loglin_impl::from_function(func_t func, range_t range_x, 
                                  size_t npoints)
-> interpol_loglin_impl
{

  auto gunc = [&func] (real_t z) {
    return func(z2x(z));
  };

  auto rgz{ rgx2rgz(range_x) };
  auto yz{
    interpol_reglin_impl::from_function(gunc, rgz, npoints)
  };
  
  return interpol_loglin_impl{std::move(yz)};
}


auto interpol_loglin_impl::transformed(func_t func) const
-> interpol_loglin_impl
{
  assert_valid();
  return interpol_loglin_impl{yz.transformed(func)};
}

auto interpol_loglin_impl::make_transform(func_t func) const
->std::shared_ptr<interpolator_impl>
{ 
  return std::make_shared<interpol_loglin_impl>(transformed(func));
}

auto interpol_loglin_impl::rescale_x(real_t scale) 
const -> interpol_loglin_impl
{       
  return interpol_loglin_impl{ yz.shift_x(log(scale)) };
}

auto interpol_loglin_impl::make_rescale_x(real_t scale) const
->std::shared_ptr<interpolator_impl>
{
  return std::make_shared<interpol_loglin_impl>(rescale_x(scale));  
}



void interpol_loglin_impl::save(datasink s) const
{
  assert_valid();
  
  s["interpolator_type"] = datastore_id;
  s["linear_interp"] = yz;
}


auto interpol_loglin_impl::from_datasource(datasource s) 
-> interpol_loglin_impl
{
  std::string styp = s["interpolator_type"];
  if (styp != datastore_id) {
    throw std::runtime_error("unexpected interpolator type in "
                             "datasource encountered");
  }
  auto yz = interpol_reglin_impl::from_datasource(s / "linear_interp");
  
  return interpol_loglin_impl{std::move(yz)};
}


} // namespace detail



auto operator*(detail::interpol_reglin_impl i, real_t a)
->detail::interpol_reglin_impl
{
  return i.transformed([a](real_t y) {return a*y;});
}


auto operator*(real_t a, detail::interpol_reglin_impl i)
->detail::interpol_reglin_impl
{
  return i * a;
}

auto operator/(detail::interpol_reglin_impl i, real_t a)
->detail::interpol_reglin_impl
{
  return i * (1.0/a);
}

auto operator/(real_t a, detail::interpol_reglin_impl i)
->detail::interpol_reglin_impl
{
  return i.transformed([a](real_t y) {return a/y;});
}


auto make_interpol_reglin(detail::interpol_reglin_impl loc)
-> interpolator 
{
  return interpolator{
      std::make_shared<detail::interpol_reglin_impl>(std::move(loc))
  };  
}

auto make_interpol_reglin(std::vector<real_t> v, interval<real_t> r)
-> interpolator 
{
  return make_interpol_reglin(
          detail::interpol_reglin_impl::from_vector(std::move(v), r));
}

auto make_interpol_reglin(std::function<real_t(real_t)> f, 
                          interval<real_t> r, size_t n)
-> interpolator 
{
  return make_interpol_reglin(
          detail::interpol_reglin_impl::from_function(f, r, n));  
}

auto make_interpol_reglin(datasource s)
-> interpolator 
{
  return make_interpol_reglin(
          detail::interpol_reglin_impl::from_datasource(s));  
}

auto make_interpol_loglin(detail::interpol_loglin_impl loc)
-> interpolator 
{
  return interpolator{
      std::make_shared<detail::interpol_loglin_impl>(std::move(loc))
  };  
}

auto make_interpol_loglin(std::vector<real_t> v, interval<real_t> r)
-> interpolator 
{
  return make_interpol_loglin(
       detail::interpol_loglin_impl::from_vector(std::move(v), r));
}

auto make_interpol_loglin(std::function<real_t(real_t)> f, 
                          interval<real_t> r, size_t n)
-> interpolator 
{
  return make_interpol_loglin(
          detail::interpol_loglin_impl::from_function(f, r, n));  
}

auto make_interpol_loglin(datasource s)
-> interpolator 
{
  return make_interpol_loglin(
          detail::interpol_loglin_impl::from_datasource(s));  
}



auto operator*(detail::interpol_loglin_impl i, real_t a)
->detail::interpol_loglin_impl
{
  return i.transformed([a](real_t y) {return a*y;});
}


auto operator*(real_t a, detail::interpol_loglin_impl i)
->detail::interpol_loglin_impl
{
  return i * a;
}

auto operator/(detail::interpol_loglin_impl i, real_t a)
->detail::interpol_loglin_impl
{
  return i * (1.0/a);
}

auto operator/(real_t a, detail::interpol_loglin_impl i)
->detail::interpol_loglin_impl
{
  return i.transformed([a](real_t y) {return a/y;});
}


} // namespace EOS_Toolkit





