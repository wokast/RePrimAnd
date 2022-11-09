#include "interpol_logspl.h"
#include <cmath>
#include <algorithm>
#include <iterator>
#include <stdexcept>

namespace EOS_Toolkit {
  
namespace detail {


const std::string interpol_logspl_impl::datastore_id {
  "cubic_monotone_spline_regular_spaced_logx"
};

auto interpol_logspl_impl::x2z(real_t x)
-> real_t
{
  return log(x);
}

auto interpol_logspl_impl::z2x(real_t z)
-> real_t
{
  return exp(z);
}

auto interpol_logspl_impl::rgz2rgx(range_t rgz)
-> range_t
{
  return {z2x(rgz.min()), z2x(rgz.max())};
  
}

auto interpol_logspl_impl::rgx2rgz(range_t rgx)
-> range_t
{
  if (rgx.min() <= 0) 
  {
    throw std::range_error("Invalid x-range for log-spaced "
                           "interpolation");
  }
  return {x2z(rgx.min()), x2z(rgx.max())};
}


interpol_logspl_impl::interpol_logspl_impl(interpol_regspl_impl yz_)
: yz{std::move(yz_)}, rgx{rgz2rgx(yz.range_x())}
{}


void interpol_logspl_impl::swap(interpol_logspl_impl& other)
{
  using std::swap;
  swap(yz, other.yz);
  swap(rgx, other.rgx);
}

void swap(interpol_logspl_impl& a, interpol_logspl_impl& b)
{ 
  a.swap(b);
}


void interpol_logspl_impl::assert_valid() const
{
  yz.assert_valid();
}

auto interpol_logspl_impl::range_x() const -> const range_t& 
{
  assert_valid();
  return rgx;
}

auto interpol_logspl_impl::range_y() const -> const range_t&
{
  assert_valid();
  return yz.range_y();
}


auto interpol_logspl_impl::operator()(real_t x) const -> real_t
{
  return yz(x2z(x));
}

auto interpol_logspl_impl::from_vector(std::vector<real_t> values, 
                                       range_t range_x) 
  -> interpol_logspl_impl
{
  auto rgz{ rgx2rgz(range_x) };
  auto yz{ interpol_regspl_impl::from_vector(std::move(values), rgz) };
  return interpol_logspl_impl{std::move(yz)};
}

auto interpol_logspl_impl::from_function(func_t func, range_t range_x, 
                                         size_t npoints)
-> interpol_logspl_impl
{

  auto gunc = [&func] (real_t z) {
    return func(z2x(z));
  };

  auto rgz{ rgx2rgz(range_x) };
       
  auto yz{
    interpol_regspl_impl::from_function(gunc, rgz, npoints)
  };
  
  return interpol_logspl_impl{std::move(yz)};
}


auto interpol_logspl_impl::transformed(func_t func) const
-> interpol_logspl_impl
{
  assert_valid();
  return interpol_logspl_impl{yz.transformed(func)};
}


auto interpol_logspl_impl::rescale_x(real_t scale) 
const -> interpol_logspl_impl
{       
  return interpol_logspl_impl{ yz.shift_x(log(scale)) };
}


auto interpol_logspl_impl::make_transform(func_t func) const
->std::shared_ptr<interpolator_impl>
{ 
  return std::make_shared<interpol_logspl_impl>(transformed(func));
}

auto interpol_logspl_impl::make_rescale_x(real_t scale) const
->std::shared_ptr<interpolator_impl>
{
  return std::make_shared<interpol_logspl_impl>(rescale_x(scale));  
}



void interpol_logspl_impl::save(datasink s) const
{
  assert_valid();
  
  s["interpolator_type"] = datastore_id;
  s["regular_spline"] = yz;
}


auto interpol_logspl_impl::from_datasource(datasource s) 
-> interpol_logspl_impl
{
  std::string styp = s["interpolator_type"];
  if (styp != datastore_id) {
    throw std::runtime_error("unexpected interpolator type in "
                             "datasource encountered");
  }
  auto yz = interpol_regspl_impl::from_datasource(
                            s / "regular_spline");
  
  return interpol_logspl_impl{std::move(yz)};
}





const std::string interpol_llogspl_impl::datastore_id {
  "cubic_monotone_spline_regular_spaced_logxlogy"
};




interpol_llogspl_impl::interpol_llogspl_impl(interpol_logspl_impl yz_)
: yz{std::move(yz_)}, 
  rgy{interpol_logspl_impl::rgz2rgx(yz.range_y())}
{}


void interpol_llogspl_impl::swap(interpol_llogspl_impl& other)
{
  using std::swap;
  swap(yz, other.yz);
  swap(rgy, other.rgy);
}

void swap(interpol_llogspl_impl& a, interpol_llogspl_impl& b)
{ 
  a.swap(b);
}


void interpol_llogspl_impl::assert_valid() const
{
  yz.assert_valid();
}

auto interpol_llogspl_impl::range_x() const -> const range_t& 
{
  assert_valid();
  return yz.range_x();
}

auto interpol_llogspl_impl::range_y() const -> const range_t&
{
  assert_valid();
  return rgy;
}


auto interpol_llogspl_impl::operator()(real_t x) const -> real_t
{
  return interpol_logspl_impl::z2x(yz(x));
}

auto interpol_llogspl_impl::from_vector(std::vector<real_t> values, 
                                        range_t range_x) 
  -> interpol_llogspl_impl
{
  std::vector<real_t> yt;
  std::transform(values.begin(), values.end(), 
                 std::back_inserter(yt), interpol_logspl_impl::x2z);
  auto yz {
    interpol_logspl_impl::from_vector(std::move(yt), range_x)
  };
  
  return interpol_llogspl_impl{std::move(yz)};
}

auto interpol_llogspl_impl::from_function(func_t func, range_t range_x, 
                                          size_t npoints)
-> interpol_llogspl_impl
{

  auto gunc = [&] (real_t z) {
    return interpol_logspl_impl::x2z(func(z));
  };
       
  auto yz{
    interpol_logspl_impl::from_function(gunc, range_x, npoints)
  };
  
  return interpol_llogspl_impl{std::move(yz)};
}


auto interpol_llogspl_impl::transformed(func_t func) const
-> interpol_llogspl_impl
{
  assert_valid();
  auto gunc { 
    [&] (real_t z) {
      return interpol_logspl_impl::x2z(
               func(interpol_logspl_impl::z2x(z))
             );
    }
  };
  return interpol_llogspl_impl{yz.transformed(gunc)};
}


auto interpol_llogspl_impl::rescale_x(real_t scale) 
const -> interpol_llogspl_impl
{       
  return interpol_llogspl_impl{ yz.rescale_x(scale) };
}


auto interpol_llogspl_impl::make_transform(func_t func) const
->std::shared_ptr<interpolator_impl>
{ 
  return std::make_shared<interpol_llogspl_impl>(transformed(func));
}

auto interpol_llogspl_impl::make_rescale_x(real_t scale) const
->std::shared_ptr<interpolator_impl>
{
  return std::make_shared<interpol_llogspl_impl>(rescale_x(scale));  
}


void interpol_llogspl_impl::save(datasink s) const
{
  assert_valid();
  
  s["interpolator_type"] = datastore_id;
  s["log_spline"] = yz;
}


auto interpol_llogspl_impl::from_datasource(datasource s) 
-> interpol_llogspl_impl
{
  std::string styp = s["interpolator_type"];
  if (styp != datastore_id) {
    throw std::runtime_error("unexpected interpolator type in "
                             "datasource encountered");
  }
  auto yz = interpol_logspl_impl::from_datasource(
                            s / "log_spline");
  
  return interpol_llogspl_impl{std::move(yz)};
}

} // namespace detail





auto make_interpol_logspl(detail::interpol_logspl_impl loc)
-> interpolator 
{
  return interpolator{
      std::make_shared<detail::interpol_logspl_impl>(std::move(loc))
  };  
}

auto make_interpol_logspl(std::vector<real_t> v, interval<real_t> r)
-> interpolator 
{
  return make_interpol_logspl(
       detail::interpol_logspl_impl::from_vector(std::move(v), r));
}

auto make_interpol_logspl(std::function<real_t(real_t)> f, 
                          interval<real_t> r, size_t n)
-> interpolator 
{
  return make_interpol_logspl(
          detail::interpol_logspl_impl::from_function(f, r, n));  
}

auto make_interpol_logspl(datasource s)
-> interpolator 
{
  return make_interpol_logspl(
          detail::interpol_logspl_impl::from_datasource(s));  
}



auto operator*(detail::interpol_logspl_impl i, real_t a)
->detail::interpol_logspl_impl
{
  return i.transformed([a](real_t y) {return a*y;});
}


auto operator*(real_t a, detail::interpol_logspl_impl i)
->detail::interpol_logspl_impl
{
  return i * a;
}

auto operator/(detail::interpol_logspl_impl i, real_t a)
->detail::interpol_logspl_impl
{
  return i * (1.0/a);
}

auto operator/(real_t a, detail::interpol_logspl_impl i)
->detail::interpol_logspl_impl
{
  return i.transformed([a](real_t y) {return a/y;});
}








auto make_interpol_llogspl(detail::interpol_llogspl_impl loc)
-> interpolator 
{
  return interpolator{
      std::make_shared<detail::interpol_llogspl_impl>(std::move(loc))
  };  
}

auto make_interpol_llogspl(std::vector<real_t> v, interval<real_t> r)
-> interpolator 
{
  return make_interpol_llogspl(
       detail::interpol_llogspl_impl::from_vector(std::move(v), r));
}

auto make_interpol_llogspl(std::function<real_t(real_t)> f, 
                           interval<real_t> r, size_t n)
-> interpolator 
{
  return make_interpol_llogspl(
        detail::interpol_llogspl_impl::from_function(f, r, n));  
}

auto make_interpol_llogspl(datasource s)
-> interpolator 
{
  return make_interpol_llogspl(
          detail::interpol_llogspl_impl::from_datasource(s));  
}



auto operator*(detail::interpol_llogspl_impl i, real_t a)
->detail::interpol_llogspl_impl
{
  return i.transformed([a](real_t y) {return a*y;});
}


auto operator*(real_t a, detail::interpol_llogspl_impl i)
->detail::interpol_llogspl_impl
{
  return i * a;
}

auto operator/(detail::interpol_llogspl_impl i, real_t a)
->detail::interpol_llogspl_impl
{
  return i * (1.0/a);
}

auto operator/(real_t a, detail::interpol_llogspl_impl i)
->detail::interpol_llogspl_impl
{
  return i.transformed([a](real_t y) {return a/y;});
}


} // namespace EOS_Toolkit





