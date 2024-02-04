#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <vector>
#include <cassert>

namespace EOS_Toolkit {



template<class C>
auto integ_forw3(const C x1, const C x2, const C x3,
                 const C y1, const C y2, const C y3) -> C
{
  C h12{ x1-x2 };
  C h23{ x2-x3 };
  C h13{ x1-x3 };
  assert((h12 < 0) && (h23 < 0));
  
  C p0{ y2 };
  C p1{ ((h23 / h12) * (y1-y2) + (h12 / h23) * (y2-y3)) / h13 };
  C p2{ (2./h13) * ((y1-y2) / h12 - (y2-y3) / h23) };
  
  return -h12 * (p0 + h12 * (p1 /2. + p2*h12/6.));
}

template<class C>
auto integ_backw3(const C x1, const C x2, const C x3,
                 const C y1, const C y2, const C y3) -> C
{
  C h12{ x1-x2 };
  C h23{ x2-x3 };
  C h13{ x1-x3 };
  assert((h12 < 0) && (h23 < 0));
  
  C p0{ y2 };
  C p1{ ((h23 / h12) * (y1-y2) + (h12 / h23) * (y2-y3)) / h13 };
  C p2{ (2./h13) * ((y1-y2) / h12 - (y2-y3) / h23) };
  
  return -h23 * (p0 + h23 * (-p1 /2. + p2*h23/6.));
}

template<class C>
auto integrate_order3(const std::vector<C>& x, const std::vector<C>& y,
                      C integ_const=0) -> std::vector<C>
{
  assert(x.size() == y.size());
  
  std::vector<C> z(x.size());
  z[0] = integ_const;
  {
    std::size_t k{ 0 };
    while (k + 2 < x.size()) {
      z[k+1] = z[k] + integ_forw3(x[k], x[k+1], x[k+2], 
                                  y[k], y[k+1], y[k+2]);
      ++k;
    }
    z[k+1] = z[k] + integ_backw3(x[k-1], x[k], x[k+1], 
                                  y[k-1], y[k], y[k+1]);
  }
  return z;
}



}

#endif
