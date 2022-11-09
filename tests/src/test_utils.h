#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include <stdexcept>
#include <string>
#include <cmath>


template<class T>
struct mapping_linear {
  typedef T vtype;
  static constexpr T y(T x) { return x;}
  static constexpr T x(T y) { return y;}
};

template<class T>
struct mapping_log {
  typedef T vtype;
  static constexpr T y(T x) { return log(x);}
  static constexpr T x(T y) { return exp(y);}
};

template<class M, class T=typename M::vtype>
class mapped_spacing {
  public:
  
  class rangeit {
    size_t i;
    const T y0,dy,x1;
    
    public:
    constexpr rangeit(size_t i_, T y0_=0., T dy_=0., T x1_=0.)
    : i{i_}, y0{y0_}, dy{dy_}, x1{x1_} {}

    T operator*() const { return std::min(x1, M::x(y0 + i*dy)); }

    rangeit& operator++() {
        ++i;
        return *this;
    }

    bool operator!=(const rangeit &o) const { return i != o.i; }    
  };
  
  constexpr mapped_spacing(T x0, T x1, size_t size)
  : itbegin{0, M::y(x0), (M::y(x1)-M::y(x0)) / size, x1}, 
    itend{size+1} {}

  rangeit begin() const { return itbegin; }
  rangeit end()   const { return itend; }
    
  private:
  rangeit itbegin,itend;
};

template<class T>
constexpr mapped_spacing< mapping_linear<T> > 
linear_spacing(T x0, T x1, size_t size) {
  return mapped_spacing< mapping_linear<T> >(x0,x1,size);
}

template<class T>
constexpr mapped_spacing< mapping_log<T> > 
log_spacing(T x0, T x1, size_t size) {
  return mapped_spacing< mapping_log<T> >(x0,x1,size);
}


class failcount {
  int count;
  std::string descr;
  public:
  failcount(std::string descr_) : count(0), descr(descr_) {}
  ~failcount();
  operator bool() const {return count ==0;}
  bool istrue(bool c, std::string msg);
  bool operator()(bool c, std::string msg) {return istrue(c,msg);}
  void fail(std::string msg);
  bool isnan(double v, std::string msg);
  bool isfinite(double v, std::string msg);
  bool isless(double smaller, double larger, std::string msg);
  bool isleq(double smaller, double larger, std::string msg);
  bool isclose(double a, double b, double reltol, double abstol, 
               std::string msg);
  void postmortem(std::string msg) const;
  
  template<typename F, typename ... A>
  bool nothrow(std::string msg, F&& f, A&& ...a);
  template<typename F, typename ... A>
  bool dothrow(std::string msg, F&& f, A&& ...a);
};

template<typename F, typename ... A>
bool failcount::nothrow(std::string msg, F&& f, A&& ...a)
{
  bool c{false};
  std::string err{"unknown exception"};
  try {
    f(std::forward<A>(a)...);
    c = true;
  }
  catch (const std::exception& e) {
    err = e.what();    
  }
  catch (...) {}
  
  istrue(c, msg+" should not throw exception");
  if (!c) postmortem(err);
  
  return c;
}

template<typename F, typename ... A>
bool failcount::dothrow(std::string msg, F&& f, A&& ...a)
{
  bool c{false};
  try {
    f(std::forward<A>(a)...);
  }
  catch (...) {
    c = true;
  }
  
  return istrue(c, msg+" should throw exception");
}




#endif
