#include <cstddef>
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
    std::size_t i;
    const T y0,dy;
    
    public:
    constexpr rangeit(std::size_t i_, T y0_=0., T dy_=0.)
    : i{i_}, y0{y0_}, dy{dy_} {}

    T operator*() const { return M::x(y0 + i*dy); }

    rangeit& operator++() {
        ++i;
        return *this;
    }

    bool operator!=(const rangeit &o) const { return i != o.i; }    
  };
  
  constexpr mapped_spacing(T x0, T x1, std::size_t size)
  : itbegin{0, M::y(x0), (M::y(x1)-M::y(x0)) / size}, 
    itend{size+1} {}

  rangeit begin() const { return itbegin; }
  rangeit end()   const { return itend; }
    
  private:
  rangeit itbegin,itend;
};

template<class T>
constexpr mapped_spacing< mapping_linear<T> > 
linear_spacing(T x0, T x1, std::size_t size) {
  return mapped_spacing< mapping_linear<T> >(x0,x1,size);
}

template<class T>
constexpr mapped_spacing< mapping_log<T> > 
log_spacing(T x0, T x1, std::size_t size) {
  return mapped_spacing< mapping_log<T> >(x0,x1,size);
}
