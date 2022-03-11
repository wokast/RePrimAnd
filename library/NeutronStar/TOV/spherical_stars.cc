#include <stdexcept>
#include <algorithm>
#include "spherical_stars_internals.h"
#include "tov_ode.h"
#include "tidal_deform_ode.h"
#include "solve_ode.h"


namespace {

auto rel_err(double a, double b) -> double
{
  return std::fabs(a - b) * 2.0 / (std::fabs(a) + std::fabs(b)); 
} 

}

namespace EOS_Toolkit {
  
auto make_tov_star(const eos_barotr eos, const real_t rho_center, 
                   const tov_acc_simple acc, const bool find_bulk, 
                   const bool find_tidal)
-> spherical_star
{
  tov_ode ode{rho_center, eos};
  tov_ode::observer obs{ode};
  auto surf{ integrate_ode_adapt(ode, acc.tov, acc.tov, 
                                 acc.minsteps, obs) };
  assert(obs.dnu.size()>0);

  auto prop{ ode.star(surf) };
  
  boost::optional<spherical_star_tidal> deform;
  if (eos.is_isentropic() and find_tidal) {
    deform = find_deform(eos, prop.center_gm1, obs.dnu, 
                         obs.rsqr, obs.lambda,
                         acc.deform);
  }

  
  auto prof = std::make_shared<details::tov_profile>(eos, prop, 
                obs.rsqr, obs.dnu, obs.lambda,
                obs.ebnd_by_r, obs.pvol_by_r);

  boost::optional<spherical_star_bulk> bulk;
  if (find_bulk) {
    bulk = find_bulk_props(*prof, acc.tov); 
  }

  return {prop, deform, bulk, prof};
}
  
auto get_tov_star_properties(const eos_barotr eos, 
                   const real_t rho_center, 
                   const tov_acc_simple acc, const bool find_bulk, 
                   const bool find_tidal)
-> spherical_star_properties
{
  tov_ode ode{rho_center, eos};
  tov_ode::observer obs{ode};
  auto surf{ integrate_ode_adapt(ode, acc.tov, acc.tov, 
                                 acc.minsteps, obs) };
  assert(obs.dnu.size()>0);

  auto prop{ ode.star(surf) };
  
  boost::optional<spherical_star_tidal> deform;
  if (eos.is_isentropic() and find_tidal) {
    deform = find_deform(eos, prop.center_gm1, obs.dnu, 
                         obs.rsqr, obs.lambda,
                         acc.deform);
  }

  boost::optional<spherical_star_bulk> bulk;
  if (find_bulk) {
    details::tov_profile prof{eos, prop, 
                obs.rsqr, obs.dnu, obs.lambda,
                obs.ebnd_by_r, obs.pvol_by_r};

    bulk = find_bulk_props(prof, acc.tov); 
  }

  return {eos, prop, deform, bulk};
}


auto make_tov_star(const eos_barotr eos, const real_t rho_center, 
                   const tov_acc_precise acc, const bool find_bulk, 
                   const bool find_tidal)
-> spherical_star
{
  tov_ode ode{rho_center, eos};
  
  auto tov_tol = [&] (const spherical_star_info&a, 
                      const spherical_star_info& b) -> bool
  {    
    return ((rel_err(a.grav_mass, b.grav_mass) < acc.mass) &&
            (rel_err(a.grav_mass+a.binding_energy, 
                     b.grav_mass+b.binding_energy) < acc.mass) &&
            (rel_err(a.circ_radius, b.circ_radius) < acc.radius) &&
            (rel_err(a.proper_volume, 
                     b.proper_volume) < 3.*acc.radius) &&
            (rel_err(a.moment_inertia, 
                     b.moment_inertia) < acc.minertia) ); 
  };
  
  auto tov_solv = [&] (real_t acc) -> spherical_star_info
  {
    return ode.star(integrate_ode_adapt(ode, acc, acc));
  }; 
  
  real_t tov_acc{ acc.mass };
  ensure_global_accuracy(tov_solv, tov_tol, tov_acc, acc.acc_min, 2.);
    
  tov_ode::observer obs{ode};
  
  auto prop{ 
    ode.star(integrate_ode_adapt(ode, tov_acc, tov_acc, 
                                 acc.minsteps, obs)) 
  };
   
  boost::optional<spherical_star_tidal> deform;
  if (eos.is_isentropic() and find_tidal) {
    auto tidal_solv = [&] (real_t acc) -> spherical_star_tidal 
    {
      return find_deform(eos, prop.center_gm1, obs.dnu, 
                         obs.rsqr, obs.lambda, acc);
    };  

    auto tidal_tol = [&] (const spherical_star_tidal&a, 
                          const spherical_star_tidal& b) -> bool
    {
      return ((rel_err(a.lambda, b.lambda) < acc.deform) && 
              (rel_err(a.k2, b.k2) < acc.deform) ); 
    };
    
    real_t tidal_acc{ acc.deform };
    deform = ensure_global_accuracy(tidal_solv, tidal_tol, 
                                    tidal_acc, acc.acc_min, 2.);

  }

  auto prof = std::make_shared<details::tov_profile>(eos, prop, 
                obs.rsqr, obs.dnu, obs.lambda,
                obs.ebnd_by_r, obs.pvol_by_r);

 
  boost::optional<spherical_star_bulk> bulk;
  if (find_bulk) {
    bulk = find_bulk_props(*prof, tov_acc); 
  }  
 
  return {prop, deform, bulk, prof};
}

auto get_tov_star_properties(const eos_barotr eos, 
                   const real_t rho_center, const tov_acc_precise acc, 
                   const bool find_bulk, const bool find_tidal)
-> spherical_star_properties
{
  tov_ode ode{rho_center, eos};
  
  auto tov_tol = [&] (const spherical_star_info&a, 
                      const spherical_star_info& b) -> bool
  {    
    return ((rel_err(a.grav_mass, b.grav_mass) < acc.mass) &&
            (rel_err(a.grav_mass+a.binding_energy, 
                     b.grav_mass+b.binding_energy) < acc.mass) &&
            (rel_err(a.circ_radius, b.circ_radius) < acc.radius) &&
            (rel_err(a.proper_volume, 
                     b.proper_volume) < 3.*acc.radius) &&
            (rel_err(a.moment_inertia, 
                     b.moment_inertia) < acc.minertia) ); 
  };
  
  auto tov_solv = [&] (real_t acc) -> spherical_star_info
  {
    return ode.star(integrate_ode_adapt(ode, acc, acc));
  }; 
  
  real_t tov_acc{ acc.mass };
  ensure_global_accuracy(tov_solv, tov_tol, tov_acc, acc.acc_min, 2.);
    
  tov_ode::observer obs{ode};
  
  auto prop{ 
    ode.star(integrate_ode_adapt(ode, tov_acc, tov_acc, 
                                 acc.minsteps, obs)) 
  };
   
  boost::optional<spherical_star_tidal> deform;
  if (eos.is_isentropic() and find_tidal) {
    auto tidal_solv = [&] (real_t acc) -> spherical_star_tidal 
    {
      return find_deform(eos, prop.center_gm1, obs.dnu, 
                         obs.rsqr, obs.lambda, acc);
    };  

    auto tidal_tol = [&] (const spherical_star_tidal&a, 
                          const spherical_star_tidal& b) -> bool
    {
      return ((rel_err(a.lambda, b.lambda) < acc.deform) && 
              (rel_err(a.k2, b.k2) < acc.deform) ); 
    };
    
    real_t tidal_acc{ acc.deform };
    deform = ensure_global_accuracy(tidal_solv, tidal_tol, 
                                    tidal_acc, acc.acc_min, 2.);

  }

  boost::optional<spherical_star_bulk> bulk;
  if (find_bulk) {
    details::tov_profile prof{eos, prop, 
                obs.rsqr, obs.dnu, obs.lambda,
                obs.ebnd_by_r, obs.pvol_by_r};

    bulk = find_bulk_props(prof, tov_acc); 
  }  
 
  return {eos, prop, deform, bulk};
}


                 
spherical_star::spherical_star(spherical_star_info info_,
                  deform_t deform_, bulk_t bulk_, pprof_t pprof_)
: spherical_star_properties(pprof_->eos(), info_, deform_, bulk_), 
  pprof{pprof_}
{
  assert(pprof);
}



auto spherical_star::profile() const 
-> const spherical_star_profile&
{
  assert(pprof);
  return *pprof;
}


auto spherical_star::nu_from_rc(real_t rc) const -> real_t
{
  return profile().nu_from_rc(rc);
}

auto spherical_star::lambda_from_rc(real_t rc) const -> real_t
{
  return profile().lambda_from_rc(rc);
}

auto spherical_star::mbary_from_rc(real_t rc) const -> real_t
{
  return profile().mbary_from_rc(rc);
}

auto spherical_star::pvol_from_rc(real_t rc) const -> real_t
{
  return profile().pvol_from_rc(rc);
}

auto spherical_star::state_from_rc(real_t rc) const 
-> eos_barotr::state
{
  return profile().state_from_rc(rc);
}

auto spherical_star::gm1_from_rc(real_t rc) const -> real_t
{
  return state_from_rc(rc).gm1();
}

auto spherical_star::rho_from_rc(real_t rc) const -> real_t
{
  return state_from_rc(rc).rho();
}

auto spherical_star::press_from_rc(real_t rc) const -> real_t
{
  return state_from_rc(rc).press();
}

auto spherical_star::eps_from_rc(real_t rc) const -> real_t
{
  return state_from_rc(rc).eps();
}

auto spherical_star::csnd_from_rc(real_t rc) const -> real_t
{
  return state_from_rc(rc).csnd();
}

auto spherical_star::ye_from_rc(real_t rc) const -> real_t
{
  return state_from_rc(rc).ye();
}

auto spherical_star::temp_from_rc(real_t rc) const -> real_t
{
  return state_from_rc(rc).temp();
}


spherical_star_properties::spherical_star_properties(
                       eos_barotr eos_, spherical_star_info info_,
                       deform_t deform_, bulk_t bulk_)
: _eos{eos_}, _info{info_}, _deform{deform_}, _bulk{bulk_} 
{}
      

auto spherical_star_properties::center_state() const 
-> eos_barotr::state
{
  return _eos.at_gm1(_info.center_gm1);
}
  
auto spherical_star_properties::center_rho() const -> real_t
{
  return center_state().rho();
}

auto spherical_star_properties::center_gm1() const -> real_t
{
  return center_state().gm1();
}

auto spherical_star_properties::center_eps() const -> real_t
{
  return center_state().eps();
}

auto spherical_star_properties::center_press() const -> real_t
{
  return center_state().press();
}

auto spherical_star_properties::center_csnd() const -> real_t
{
  return center_state().csnd();
}

auto spherical_star_properties::center_ye() const -> real_t
{
  return center_state().ye();
}


auto spherical_star_properties::grav_mass() const -> real_t
{
  return _info.grav_mass;
}

auto spherical_star_properties::bary_mass() const -> real_t
{
  return _info.binding_energy + _info.grav_mass;
}

auto spherical_star_properties::binding_energy() const -> real_t
{
  return _info.binding_energy;
}

auto spherical_star_properties::circ_radius() const -> real_t
{
  return _info.circ_radius;
}

auto spherical_star_properties::proper_volume() const -> real_t
{
  return _info.proper_volume;
}

auto spherical_star_properties::moment_inertia() const -> real_t
{
  return _info.moment_inertia;
}

auto spherical_star_properties::has_deform() const -> bool
{
  return bool(_deform);
}

auto spherical_star_properties::deformability() const 
-> const spherical_star_tidal&
{
  return _deform.value();
}

auto spherical_star_properties::has_bulk() const -> bool
{
  return bool(_bulk);
}

auto spherical_star_properties::bulk() const 
-> const spherical_star_bulk&
{
  return _bulk.value();
}

auto spherical_star_properties::eos() const -> const eos_barotr&
{
  return _eos;
}



} // namespace EOS_Toolkit



 
