#include "interpol.h"
#include <gsl/gsl_spline.h>

namespace EOS_Toolkit {
namespace detail {

struct wrap_interp_accel {
  gsl_interp_accel *p;
  wrap_interp_accel();
  ~wrap_interp_accel();
};

struct wrap_interp_cspline {
  gsl_interp *p;
  wrap_interp_cspline(const std::vector<double>& x, 
                     const std::vector<double>& y);
  ~wrap_interp_cspline();
};


///Implementation belonging to interface class \ref cspline_mono
class cspline_mono_impl {  
  public:
  using range_t = cspline_mono::range_t;
  
  cspline_mono_impl()            = delete;
  cspline_mono_impl(const cspline_mono_impl&) = delete;
  cspline_mono_impl(cspline_mono_impl&&)      = delete;
  cspline_mono_impl& operator=(cspline_mono_impl&&)      = delete;
  cspline_mono_impl& operator=(const cspline_mono_impl&) = delete;
  
  cspline_mono_impl(std::vector<double> x_, std::vector<double> y_);
       
  const range_t& range_x() const {return rgx;}
  const range_t& range_y() const {return rgy;}
  real_t operator()(real_t t) const;
  
  private:
  std::vector<double> x,y; 
  //order matters, x,y, must come before interp 
  mutable wrap_interp_accel acc;
  wrap_interp_cspline interp;
  
  range_t rgx;
  range_t rgy;
};

}
}

