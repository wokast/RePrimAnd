#include "interpol_regspl.h"
#include <cmath>
#include <algorithm>
#include <iterator>
#include <stdexcept>

namespace EOS_Toolkit {
  
namespace detail {


const std::string interpol_regspl_impl::datastore_id {
  "cubic_monotone_spline_regular_spaced"
};

interpol_regspl_impl::interpol_regspl_impl(
                             interpol_regspl_impl&& other)
{
  other.swap(*this);
}

interpol_regspl_impl& interpol_regspl_impl::operator=(
                                        interpol_regspl_impl other)
{
  other.swap(*this);
  return *this;
}

void interpol_regspl_impl::swap(interpol_regspl_impl& other)
{
  using std::swap;
  swap(segs, other.segs);
  swap(rgx, other.rgx);
  swap(rgy, other.rgy);
  swap(dx, other.dx);
}

void swap(interpol_regspl_impl& a, interpol_regspl_impl& b)
{ 
  a.swap(b);
}

void interpol_regspl_impl::assert_valid() const
{
  assert(!segs.empty());
}


auto interpol_regspl_impl::get_dx(const range_t& rgx, 
                                     std::size_t nsegs)
-> real_t
{
  if (nsegs < 2) {
    throw std::range_error("interpol_regspl_impl: need as least 3 "
                           "sample points");
  }
  if (rgx.length() <= 0) {
    throw std::range_error("interpol_regspl_impl: degenerate x-range");
  }
  
  return rgx.length() / nsegs;
}

auto interpol_regspl_impl::get_rgy(const std::vector<real_t>& y)
-> range_t
{
  auto ext = std::minmax_element(y.begin(), y.end());
  return {*ext.first, *ext.second};
}

interpol_regspl_impl::interpol_regspl_impl(
                  std::vector<segment> segments, 
                  range_t range_x, range_t range_y)
: segs{std::move(segments)}, rgx{range_x}, rgy{range_y}, 
  dx{get_dx(rgx, segs.size())}
{}

auto interpol_regspl_impl::range_x() const -> const range_t& 
{
  assert_valid();
  return rgx;
}

auto interpol_regspl_impl::range_y() const -> const range_t&
{
  assert_valid();
  return rgy;
}


auto interpol_regspl_impl::segment::operator()(real_t t) const 
-> real_t
{
  return ((c[0] * t + c[1]) * t + c[2]) * t + c[3];
}

auto interpol_regspl_impl::segment::hermite(real_t y0, real_t y1, 
                                            real_t m0, real_t m1) 
-> segment
{
  real_t delta{ y1 - y0 };
  return segment{{
    -2. * delta + m0 + m1,
    3. * delta - 2. * m0 - m1,
    m0, 
    y0
  }};
}

auto interpol_regspl_impl::make_seg(std::array<real_t, 4> y) 
-> segment
{
  using std::max;
  using std::min;
  
  std::array<real_t, 4> delta, maxd, mind, lm;
  
  for (std::size_t i=0; i < y.size()-1; ++i) {
    delta[i] = y[i+1] - y[i];
    maxd[i] = 3.0 * max(0., delta[i]);
    mind[i] = 3.0 * min(0., delta[i]);
  }

  for (std::size_t i=1; i < y.size()-1; ++i) {
    real_t m_i{ (delta[i-1] + delta[i]) / 2. };
    m_i = max( m_i, max(mind[i-1], mind[i]) ); 
    lm[i] = min( m_i, min(maxd[i-1], maxd[i]) ); 
  }
  
  return segment::hermite(y[1], y[2], lm[1], lm[2]);
   
}

auto interpol_regspl_impl::from_vector(std::vector<real_t> y, 
                                       range_t range_x) 
-> interpol_regspl_impl
{
  std::size_t n{ y.size() };
  auto range_y{ get_rgy(y) };
  //~ real_t dx{ range_x.length() / (n-1) };
  
  std::vector<segment> segs;
  
  segs.push_back(make_seg({y[0] - (y[1]-y[0]), 
                          y[0], y[1], y[2]}));
  
  for(std::size_t i=1; i < n - 2; ++i)
  {
    segs.push_back(make_seg({y[i-1], y[i], y[i+1], y[i+2]}));
  }
  
  segs.push_back(make_seg({y[n-3], y[n-2], y[n-1], 
                          y[n-1] + (y[n-1]-y[n-2])}));
  assert(segs.size() + 1 == y.size());
   
  return interpol_regspl_impl{std::move(segs), range_x, range_y};
}
  
auto interpol_regspl_impl::from_function(func_t func, range_t range_x, 
                                         size_t npoints)
-> interpol_regspl_impl
{
  real_t dx = get_dx(range_x, npoints-1);
  std::vector<real_t> y;
  for (std::size_t k=0; k<npoints; ++k) {
    real_t x{ range_x.min() + dx*k }; 
    y.push_back(func(range_x.limit_to(x)));
  }
  return from_vector(std::move(y), range_x);
}

auto interpol_regspl_impl::transformed(func_t func) const
-> interpol_regspl_impl
{
  assert_valid();

  auto gunc = [&] (real_t x) {
    return func((*this)(x));
  };

  return from_function(gunc, rgx, segs.size()+1);
}

auto interpol_regspl_impl::rescale_x(real_t scale) const 
-> interpol_regspl_impl
{
  assert_valid();
  
  range_t rgxnew{ rgx.min() * scale, rgx.max() * scale };

  auto gunc = [&] (real_t xnew) {
    return (*this)(xnew / scale);
  };

  return from_function(gunc, rgxnew, segs.size()+1);
}

auto interpol_regspl_impl::shift_x(real_t offset) const 
-> interpol_regspl_impl
{
  assert_valid();
  
  range_t rgxnew{ rgx.min() + offset, rgx.max() + offset };

  auto gunc = [&] (real_t xnew) {
    return (*this)(xnew - offset);
  };

  return from_function(gunc, rgxnew, segs.size()+1);
}


auto interpol_regspl_impl::make_transform(func_t func) const
-> std::shared_ptr<interpolator_impl>
{ 
  return std::make_shared<interpol_regspl_impl>(transformed(func));
}

auto interpol_regspl_impl::make_rescale_x(real_t scale) const
-> std::shared_ptr<interpolator_impl>
{
  return std::make_shared<interpol_regspl_impl>(rescale_x(scale));  
}



auto interpol_regspl_impl::from_datasource(datasource s) 
-> interpol_regspl_impl
{
  std::string styp = s["interpolator_type"];
  if (styp != datastore_id) {
    throw std::runtime_error("unexpected interpolator type in "
                             "datasource encountered");
  }
  std::vector<real_t> y = s["sample_values"];
  interval<real_t> rg = s["range_x"];
 
  return from_vector(std::move(y), rg);
}

void interpol_regspl_impl::save(datasink s) const
{
  assert_valid();
  
  std::vector<real_t> y;
  for (auto s : segs) 
  {
    y.push_back(s(0.));
  }
  y.push_back(segs.back()(1.));
  
  s["interpolator_type"] = datastore_id;
  s["sample_values"] = y;
  s["range_x"] = rgx;
}


/**
If x is outside the tabulated range, the function value at the 
closest boundary is returned
*/
real_t interpol_regspl_impl::operator()(real_t x) const
{
  assert_valid();
  
  real_t i{ (x - rgx.min()) / dx };
  real_t j{ std::max(0., floor(i)) };
  std::size_t k{ std::min(std::size_t(j), segs.size()-1) };
  
  return segs[k](i - k);
}



} // namespace detail


auto make_interpol_regspl(
         detail::interpol_regspl_impl loc)
-> interpolator 
{
  return interpolator{
      std::make_shared<detail::interpol_regspl_impl>(std::move(loc))
  };  
}

auto make_interpol_regspl(std::vector<real_t> v, interval<real_t> r)
-> interpolator 
{
  return make_interpol_regspl(
          detail::interpol_regspl_impl::from_vector(std::move(v), r));
}

auto make_interpol_regspl(std::function<real_t(real_t)> f, 
                          interval<real_t> r, std::size_t n)
-> interpolator 
{
  return make_interpol_regspl(
          detail::interpol_regspl_impl::from_function(f, r, n));  
}

auto make_interpol_regspl(datasource s)
-> interpolator 
{
  return make_interpol_regspl(
          detail::interpol_regspl_impl::from_datasource(s));  
}




auto operator*(detail::interpol_regspl_impl i, real_t a)
->detail::interpol_regspl_impl
{
  return i.transformed([a](real_t y) {return a*y;});
}


auto operator*(real_t a, detail::interpol_regspl_impl i)
->detail::interpol_regspl_impl
{
  return i * a;
}

auto operator/(detail::interpol_regspl_impl i, real_t a)
->detail::interpol_regspl_impl
{
  return i * (1.0/a);
}

auto operator/(real_t a, detail::interpol_regspl_impl i)
->detail::interpol_regspl_impl
{
  return i.transformed([a](real_t y) {return a/y;});
}





} // namespace EOS_Toolkit





