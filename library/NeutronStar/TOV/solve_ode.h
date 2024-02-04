#ifndef SOLVE_ODE_H
#define SOLVE_ODE_H

#include <algorithm>
#include <functional>
#include <cmath>
#include <cassert>
#include <boost/numeric/odeint.hpp>


namespace EOS_Toolkit {



template<class ODE>
class collecting_observer {
  public:
  using state_t = typename ODE::state_t;
  using value_t = typename state_t::value_type;
  using vec_t   = std::vector<value_t>;
  using data_t  = std::array<vec_t, std::tuple_size<state_t>::value>;
  
  vec_t x;
  data_t y;
  
  collecting_observer() = default;

  void operator()(const state_t& snew, value_t xnew) 
  {
    x.push_back(xnew);
    for (std::size_t i=0; i<y.size(); ++i) y[i].push_back(snew[i]);  
  }

};


template<class ODE>
class no_observer {
  public:
  using state_t = typename ODE::state_t;
  using value_t = typename state_t::value_type;
  
  no_observer() = default;

  void operator()(const state_t& snew, value_t xnew) 
  {}
};



template<class ODE, class OBS=no_observer<ODE>, 
         class S=typename ODE::state_t, 
         class R=typename ODE::value_t>
auto integrate_ode_adaptive(const ODE& ode, const S& s0, 
                            const R x0, const R x1, const real_t acc,  
                            const std::size_t nsample, OBS&& obs=OBS()) 
-> S
{
  using stepper_t = boost::numeric::odeint::runge_kutta_cash_karp54<S>;
  
  assert(std::isfinite(x0));
  assert(std::isfinite(x1));
  assert(nsample>1);
  
  S s{ s0 };
  R dx{ (x1-x0) / nsample };
  
  boost::numeric::odeint::integrate_n_steps( 
    boost::numeric::odeint::make_controlled<stepper_t>(acc, acc),
    std::ref(ode), s , x0 , dx , nsample , std::ref(obs)
  ); 
  
   
  return s;
}

template<class ODE, class OBS=no_observer<ODE>, class S=typename ODE::state_t, 
         class R=typename ODE::value_t>
auto integrate_ode_adaptive(const ODE& ode, const real_t acc, 
                            const std::size_t nsample,
                            OBS&& obs=OBS()) 
-> S
{
  R x0{ ode.x_start() };
  R x1{ ode.x_end() };
  
  return integrate_ode_adaptive(ode, ode.initial_data(), 
                             x0, x1, acc, nsample, obs);
}


template<class ODE, class S=typename ODE::state_t, 
         class R=typename ODE::value_t>
auto integrate_ode_fixed(const ODE& ode, const S& s0, 
  const R x0, const R x1, const std::size_t nsample) 
-> S
{
  using stepper_t = boost::numeric::odeint::runge_kutta_cash_karp54<S>;

  assert(std::isfinite(x0));
  assert(std::isfinite(x1));
  assert(nsample>1);
  
  stepper_t stepper;
  S s{ s0 };
  R dx0{ (x1-x0) / nsample };
  
  boost::numeric::odeint::integrate_n_steps(
    stepper,
    std::ref(ode), s, x0, dx0, nsample);    
   
  return s;
}


template<class ODE, class OBS, class S=typename ODE::state_t, 
         class R=typename ODE::value_t>
auto integrate_ode_fixed(const ODE& ode, const S& s0, 
  const R x0, const R x1,
  const std::size_t nsample, OBS& observer) 
-> S
{
  using stepper_t = boost::numeric::odeint::runge_kutta_cash_karp54<S>;

  assert(std::isfinite(x0));
  assert(std::isfinite(x1));
  assert(nsample>1);

  stepper_t stepper;
  S s{ s0 };
  R dx0{ (x1-x0) / nsample };
  
  boost::numeric::odeint::integrate_n_steps(
    stepper,
    std::ref(ode), s, x0, dx0, nsample, std::ref(observer));    
   
  return s;
}



template<class ODE, class OBS, class S=typename ODE::state_t, 
         class R=typename ODE::value_t>
auto integrate_ode_fixed(const ODE &ode, 
                         const std::size_t nsample,
                         OBS& observer) -> S
{
  R x0{ ode.x_start() };
  R x1{ ode.x_end() };
  return integrate_ode_fixed(ode, ode.initial_data(), 
                             x0, x1,
                             nsample, observer);
}


template<class ODE, class S=typename ODE::state_t, 
         class R=typename ODE::value_t>
auto integrate_ode_fixed(const ODE &ode, 
                         const std::size_t nsample) -> S
{
  R x0{ ode.x_start() };
  R x1{ ode.x_end() };
  return integrate_ode_fixed(ode, ode.initial_data(), 
                             x0, x1,
                             nsample);
}




}

#endif
