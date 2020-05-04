#ifndef INTERPOL_H
#define INTERPOL_H
#include <functional>
#include <memory>
#include <vector>
#include "config.h"
#include "intervals.h"

namespace EOS_Toolkit {


namespace detail {
  class cspline_mono_impl;
}


///Lookup table
/**
Approximates arbitrary functions over a given range using fast linear 
interpolation. 
**/
class lookup_table {
  public:

  using func_t  = std::function<real_t(real_t)>;
  using range_t = interval<real_t>;

  ///Default constructor. 
  lookup_table()                    = default;
  ~lookup_table()                   = default;
  lookup_table(const lookup_table&) = default;
  lookup_table(lookup_table&&)      = default;
  lookup_table& operator=(const lookup_table&) = default;
  lookup_table& operator=(lookup_table&&) = default;
  ///Sample from function
  lookup_table(func_t func, range_t range, size_t npoints);

  ///Valid range. 
  const range_t &range_x() const {return rgx;}

  ///Value range
  const range_t &range_y() const {return rgy;}

  ///Look up value. 
  real_t operator()(real_t x) const;


  private:

  std::vector<real_t> y{0,0};
  real_t dxinv{0.0};
  range_t rgx{0,0};
  range_t rgy{0,0};
};



///Lookup table designed to cover many orders of magnitude in x
class lookup_table_magx{
  public:

  using func_t = lookup_table::func_t;
  using range_t = lookup_table::range_t;

  lookup_table_magx()                           = default;
  ~lookup_table_magx()                          = default;
  lookup_table_magx(const lookup_table_magx&) = default;
  lookup_table_magx(lookup_table_magx&&)      = default;
  lookup_table_magx& operator=(const lookup_table_magx&) = default;
  lookup_table_magx& operator=(lookup_table_magx&&) = default;

  
  ///Sample from function
  lookup_table_magx(func_t func, range_t range, size_t npoints, 
                    int magnitudes);
  
  ///Valid range. 
  const range_t &range_x() const {return rgx;}

  ///Value range
  const range_t &range_y() const {return tbl.range_y();}

  ///Look up value. 
  real_t operator()(real_t x) const;
  
  private:
  lookup_table tbl;  
  range_t rgx{0,0};
  real_t x_offs{1.0};
};


///A monotonic cubic spline.
class cspline_mono {
  std::shared_ptr<const detail::cspline_mono_impl> pimpl;
  const detail::cspline_mono_impl& valid() const;
  
  public:
  using range_t = interval<real_t>;
  
  cspline_mono(const std::vector<real_t>& x, 
               const std::vector<real_t>& y);
  cspline_mono()                                = default;
  cspline_mono(const cspline_mono&)             = default;
  cspline_mono(cspline_mono&&)                  = default;
  cspline_mono& operator=(const cspline_mono&)  = default;
  cspline_mono& operator=(cspline_mono&&)       = default;
  ~cspline_mono()                               = default;

  const range_t& range_x() const;
  const range_t& range_y() const;
  real_t operator()(real_t x) const;
  
};

}

#endif

