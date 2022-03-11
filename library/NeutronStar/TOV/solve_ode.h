#ifndef SOLVE_ODE_H
#define SOLVE_ODE_H

#include <algorithm>
#include <functional>
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


template<class ODE, class OBS, class S=typename ODE::state_t, 
         class R=typename ODE::value_t>
auto integrate_ode_adapt(const ODE& ode, const S& s0, 
  const R x0, const R x1,
  const R acc_abs, const R acc_rel, 
  const std::size_t nsample, OBS& observer) 
-> S
{
  using stepper_t = boost::numeric::odeint::runge_kutta_cash_karp54<S>;

  assert(nsample > 0);
  S s{ s0 };
  R dx0{ (x1-x0) / nsample };
  
  boost::numeric::odeint::integrate_const(
    boost::numeric::odeint::make_controlled<stepper_t>(acc_abs, acc_rel),
    std::ref(ode), s, x0, x1, dx0, std::ref(observer));    
   
  return s;
}

template<class ODE, class S=typename ODE::state_t, 
         class R=typename ODE::value_t>
auto integrate_ode_adapt(const ODE& ode, const S& s0, 
  const R x0, const R x1, const R dx0,
  const R acc_abs, const R acc_rel) 
-> S
{
  using stepper_t = boost::numeric::odeint::runge_kutta_cash_karp54<S>;

  S s{ s0 };        
  boost::numeric::odeint::integrate_adaptive(
    boost::numeric::odeint::make_controlled<stepper_t>(acc_abs, acc_rel),
    std::ref(ode), s, x0, x1, dx0);    

  return s;
}




template<class ODE, class OBS, class S=typename ODE::state_t, 
         class R=typename ODE::value_t>
auto integrate_ode_adapt(const ODE &ode, 
                         const R acc_abs,  const R acc_rel, 
                         const std::size_t nsample,
                         OBS& observer) -> S
{
  R x0{ ode.x_start() };
  R x1{ ode.x_end() };
  return integrate_ode_adapt(ode, ode.initial_data(), 
                             x0, x1, acc_abs, acc_rel, 
                             nsample, observer);
}

template<class ODE, class S=typename ODE::state_t, 
         class R=typename ODE::value_t>
auto integrate_ode_adapt(const ODE &ode, const R acc_abs,  
                         const R acc_rel) -> S
{
  const R x0{ ode.x_start() };
  const R x1{ ode.x_end() };
  const R dx0{ (x1-x0) * 1e-3 };
  return integrate_ode_adapt(ode, ode.initial_data(), 
                             x0, x1, dx0, acc_abs, acc_rel);
}


template<class V, class F, class E>
auto ensure_global_accuracy(const F& func, const E& tol, V& acc, 
                            const V acc_min, const V ref_fac=2.) 
-> decltype( func(acc) )
{
  auto f0{ func(acc) };
  
  for( bool cont=true; cont;) {
    acc /= ref_fac;
    if (acc<acc_min) {
      throw std::runtime_error("Could not ensure desired accuracy");
    }
    auto f1{ func(acc) };
    cont = !tol(f1,f0);
    f0 = f1;
  }
  return f0;
}


}

#endif
