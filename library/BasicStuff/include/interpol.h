#ifndef INTERPOL_H
#define INTERPOL_H
#include <functional>
#include <memory>
#include <vector>
#include "config.h"
#include "intervals.h"
#include "datastore.h"

namespace EOS_Toolkit {


namespace detail {


class interpolator_impl {
  public:

  using func_t  = std::function<real_t(real_t)>;
  
  using range_t = interval<real_t>;
  
  ///Valid range. 
  virtual const range_t &range_x() const =0;
  
  ///Value range
  virtual const range_t &range_y() const =0;

  ///Evaluate. 
  virtual real_t operator()(real_t x) const =0;

  virtual void save(datasink dst) const =0;
    
  virtual ~interpolator_impl() = default;
    
  bool contains(real_t x) const {return range_x().contains(x);}
  bool contains(range_t r) const {return range_x().contains(r);}
  
  virtual auto make_transform(func_t) const
  ->std::shared_ptr<interpolator_impl> =0;

  virtual auto make_rescale_x(real_t scale) const
  ->std::shared_ptr<interpolator_impl> =0;

};
}

class interpolator : public detail::interpolator_impl {
  std::shared_ptr<const detail::interpolator_impl> pimpl;
  const detail::interpolator_impl& valid() const;
  
  public:
  using range_t = interval<real_t>;
  using interpolator_impl::func_t;
  
  interpolator(std::shared_ptr<interpolator_impl> impl_);
  interpolator()                                = default;
  interpolator(const interpolator&)             = default;
  interpolator(interpolator&&)                  = default;
  interpolator& operator=(const interpolator&)  = default;
  interpolator& operator=(interpolator&&)       = default;
  ~interpolator()                               = default;

  auto range_x() const -> const range_t& final;
  
  auto range_y() const ->const range_t& final;
  
  auto operator()(real_t x) const -> real_t final 
  {
    return valid()(x);
  }
  
  void save(datasink s) const final {valid().save(s);}
  
  auto transformed(func_t) const ->interpolator;
  
  auto make_transform(func_t) const 
  ->std::shared_ptr<interpolator_impl> final;

  auto rescale_x(real_t scale) const -> interpolator;

  auto make_rescale_x(real_t scale) const 
  ->std::shared_ptr<interpolator_impl> final;

};

auto make_interpolator(datasource) 
->interpolator;

auto make_interpol_reglin(std::vector<real_t>, interval<real_t>) 
->interpolator;

auto make_interpol_reglin(std::function<real_t(real_t)>, 
                          interval<real_t>, size_t)
->interpolator;

auto make_interpol_reglin(datasource) 
->interpolator;

auto make_interpol_loglin(std::vector<real_t>, interval<real_t>) 
->interpolator;

auto make_interpol_loglin(std::function<real_t(real_t)>, 
                          interval<real_t>, size_t)
->interpolator;

auto make_interpol_loglin(datasource) 
->interpolator;


auto make_interpol_regspl(std::vector<real_t> v, 
                                   interval<real_t> r)
-> interpolator;

auto make_interpol_regspl(std::function<real_t(real_t)> f, 
                          interval<real_t> r, size_t n)
-> interpolator;

auto make_interpol_regspl(datasource s)
-> interpolator;


auto make_interpol_logspl(std::vector<real_t>, interval<real_t>) 
-> interpolator;

auto make_interpol_logspl(std::function<real_t(real_t)>, 
                          interval<real_t>, size_t)
-> interpolator;

auto make_interpol_logspl(datasource) 
-> interpolator;





auto make_interpol_llogspl(std::vector<real_t>, interval<real_t>) 
-> interpolator;

auto make_interpol_llogspl(std::function<real_t(real_t)>, 
                           interval<real_t>, size_t)
-> interpolator;

auto make_interpol_llogspl(datasource) 
-> interpolator;




auto make_interpol_pchip_spline(std::vector<real_t> x, 
                                std::vector<real_t> y)
-> interpolator;

auto make_interpol_pchip_spline(std::vector<real_t> x, 
                                std::function<real_t(real_t)> f)
-> interpolator;

auto make_interpol_pchip_spline(datasource s)
-> interpolator;



auto operator*(real_t, interpolator)
->interpolator;

auto operator*(interpolator, real_t)
->interpolator;

auto operator/(interpolator, real_t)
->interpolator;

auto operator/(real_t, interpolator)
->interpolator;


namespace detail {
  

template<> struct source_proxy_reader<interpolator> {
  static void read(const datasource& s, std::string n, 
                   interpolator& t) 
  {
    t = make_interpolator(s / n);
  }
};

  
template<> struct source_proxy_writer<interpolator> {
  static void write(datasink& s, std::string n, 
                    const interpolator& t) 
  {
    t.save(s / n);
  }
};
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




}

#endif

