#ifndef INTERVALS_H
#define INTERVALS_H

#include <cmath>
#include <algorithm>
#include <cassert>

namespace EOS_Toolkit {

/**Class representing a closed interval

@tparam T the underlying numerical type, e.g. double
**/
template<class T>
class interval {
  T min_{0};  ///< Lower boundary 
  T max_{0};  ///< Upper boundary
  public:
  
  ///Default constructor yields empty range \f$ [0,0] \f$.
  interval() noexcept = default;
  
  ~interval() noexcept = default;
  interval(const interval&) noexcept = default;
  interval(interval&&) noexcept = default;
  interval& operator=(const interval&) noexcept = default;
  interval& operator=(interval&&) noexcept = default;
  
  /**\brief Construct from minimum and maximum
  
  @param min Lower boundary of interval
  @param max Upper boundary of interval
  
  \pre Aborts unless max >= min 
  \pre Aborts unless max and min are both finite values
  **/
  interval(T min, T max) : min_(min), max_(max) 
  {
    assert(min_ <= max_);
    assert(std::isfinite(min_));
    assert(std::isfinite(max_));
  }
  
  /**
  @param x Value to test
  @return If value is contained in closed interval \f$[a,b]\f$. 
  
  \note Returns false if x is NAN or INF.
  **/
  bool contains(const T& x) const 
  {
    return (x >= min_) && (x <= max_);
  }
  
  /**\brief Limit value to interval
  
  @param x Value to constrain
  @return x if it inside interval, else the closest boundary
  **/
  T limit_to(const T& x) const 
  {
    return std::min(std::max(min_, x), max_);
  }
  
  /**@return Lower boundary **/
  T min() const noexcept {return min_;}
  
  /**@return Upper boundary **/
  T max() const noexcept {return max_;}
  
  /** @return Length of interval **/
  T length() const {return max_ - min_;}
};

///Test if value is above interval
template<class T>
bool operator>(T x, const interval<T>& i) 
{
  return x > i.max();  
}

///Test if value is below interval
template<class T>
bool operator<(T x, const interval<T>& i) 
{
  return x < i.min();
}



}// namespace EOS_Toolkit 

#endif
