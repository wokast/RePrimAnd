#include <stdexcept>
#include <algorithm>
#include "spherical_stars_internals.h"
#include "tov_ode.h"
#include "tidal_deform_ode.h"
#include "solve_ode.h"



namespace EOS_Toolkit {
  
  
  
star_accuracy_spec::star_accuracy_spec(real_t acc_mass_, 
                   real_t acc_radius_, real_t acc_minertia_, 
                   std::size_t minsteps_, bool need_deform_,
                   real_t acc_deform_, bool need_bulk_)
: acc_mass{acc_mass_}, acc_radius{acc_radius_}, 
  acc_minertia{acc_minertia_}, minsteps{minsteps_}, 
  need_deform{need_deform_}, acc_deform{acc_deform_},
  need_bulk{need_bulk_} 
  {
  if (!(acc_mass > 0)) {
    throw std::runtime_error("Tolerance for mass must be greater zero");  
  }
  if (!(acc_radius > 0)) {
    throw std::runtime_error("Tolerance for radius must be "
                             "greater zero");  
  }
  if (!(acc_minertia > 0)) {
    throw std::runtime_error("Tolerance for moment of inertia "
                             " must be greater zero");  
  }
  if (!(acc_deform > 0)) {
    throw std::runtime_error("Tolerance for tidal deformability "
                             " must be greater zero");  
  }
}  

auto star_acc_simple(bool need_deform, 
                     bool need_bulk, 
                     real_t acc_tov,
                     real_t acc_deform, 
                     std::size_t minsteps) 
-> star_accuracy_spec
{
  return star_accuracy_spec(acc_tov, acc_tov, acc_tov, minsteps,
                            need_deform, acc_deform, need_bulk);
}

auto star_acc_detailed(bool need_deform,
                       bool need_bulk, 
                       real_t acc_mass, 
                       real_t acc_radius, 
                       real_t acc_minertia,
                       real_t acc_deform, 
                       std::size_t minsteps) 
-> star_accuracy_spec
{
  return star_accuracy_spec(acc_mass, acc_radius, acc_minertia,
                    minsteps, need_deform, acc_deform, need_bulk);
}


auto tov_solver_adaptive::get_star_properties(const eos_barotr& eos, 
       real_t rho_center, const star_accuracy_spec acc) 
-> spherical_star_properties
{
  return engine::get_star_properties(eos, rho_center, 
                                     heuristic_params_accuracy(acc));  
} 

auto tov_solver_adaptive::get_star(const eos_barotr& eos, 
       real_t rho_center, const star_accuracy_spec acc)
-> spherical_star
{
  return engine::get_star(eos, rho_center, 
                          heuristic_params_accuracy(acc));  
}

auto tov_solver_adaptive::heuristic_params_accuracy(
                                const star_accuracy_spec& acc)
-> engine::parameters
{
  if (acc.acc_mass < 1e-9) 
  {
    throw std::runtime_error("Requested accuracy for mass too high. "
    "Adaptive TOV solver only calibrated for tolerances >= 1e-9");
  }
  if (acc.acc_radius < 1e-9) 
  {
    throw std::runtime_error("Requested accuracy for radius too high. "
    "Adaptive TOV solver only calibrated for tolerances >= 1e-9");
  }
  if (acc.acc_minertia < 1e-9) 
  {
    throw std::runtime_error("Requested accuracy for moment of inertia " 
      "too high. Adaptive TOV solver only calibrated for "
      "tolerances >= 1e-9");
  }
  if (acc.acc_deform < 1e-9) 
  {
    throw std::runtime_error("Requested accuracy for tidal deformability " 
      "too high. Adaptive TOV solver only calibrated for "
      "tolerances >= 1e-9");
  }
  
  const real_t wdiv_tidal_default { 0.9 };
  std::size_t nsamp_tov{ acc.minsteps };
  if (acc.need_deform) 
  {
    const real_t n{ 
      std::max(50., 1e4 / pow(acc.acc_deform/1e-6, 0.5))
    };
    if (nsamp_tov < n) nsamp_tov = n;
  }
  const real_t acc_tov { 
    std::min({acc.acc_mass/100., acc.acc_radius/20., 
              acc.acc_minertia / 200.})
  };
  const std::size_t nsamp_tidal { 
    std::max(std::size_t{20}, nsamp_tov/2) 
  };
  const real_t acc_tidal { acc.acc_deform / 1e2 };
  
  return {acc.need_bulk, acc.need_deform, nsamp_tov, acc_tov,
          nsamp_tidal, acc_tidal, wdiv_tidal_default, acc.acc_radius};
}

auto tov_solver_fixstep::get_star_properties(const eos_barotr& eos, 
       real_t rho_center, const star_accuracy_spec acc) 
-> spherical_star_properties
{
  return engine::get_star_properties(eos, rho_center, 
                                     heuristic_params_accuracy(acc));  
} 

auto tov_solver_fixstep::get_star(const eos_barotr& eos, 
       real_t rho_center, const star_accuracy_spec acc)
-> spherical_star
{
  return engine::get_star(eos, rho_center, 
                          heuristic_params_accuracy(acc));  
}

auto tov_solver_fixstep::heuristic_params_accuracy(
                                const star_accuracy_spec& acc)
-> engine::parameters
{
  if (acc.acc_mass < 1e-9) 
  {
    throw std::runtime_error("Requested accuracy for mass too high. "
    "Fix-step TOV solver only calibrated for tolerances >= 1e-9");
  }
  if (acc.acc_radius < 1e-9) 
  {
    throw std::runtime_error("Requested accuracy for radius too high. "
    "Fix-step TOV solver only calibrated for tolerances >= 1e-9");
  }
  if (acc.acc_minertia < 1e-9) 
  {
    throw std::runtime_error("Requested accuracy for moment of inertia " 
      "too high. Fix-step TOV solver only calibrated for "
      "tolerances >= 1e-9");
  }
  if (acc.acc_deform < 1e-9) 
  {
    throw std::runtime_error("Requested accuracy for tidal deformability " 
      "too high. Fix-step TOV solver only calibrated for "
      "tolerances >= 1e-9");
  }
  
  
  
  const std::size_t nsub_tidal_default{ 2 };
  const real_t wdiv_tidal_default { 0.9 };
  
  const std::array<real_t,6> order {{1.8, 1.8, 1.8, 1.8, 1.8, 1.6}};
  const std::array<real_t,6> resr {{680, 700, 190, 380, 600, 2600}};
  const real_t tolr{ 1e-5 };
  const std::array<real_t,6> tol {{
    acc.acc_mass, acc.acc_mass, 
    acc.acc_radius, 3*acc.acc_radius, //volume = length^3
    acc.acc_minertia, acc.acc_deform
  }};
  std::array<bool,6> needed{
    true, true, true, true, true, 
    acc.need_deform
  };
  real_t mres{ static_cast<real_t>(acc.minsteps) };
  for (std::size_t i=0; i<order.size(); ++i) 
  {
    assert(tol[i] > 0);
    assert(order[i] > 0);
    if (needed[i]) 
    {
      real_t ires{ resr[i] * pow( tolr / tol[i], 1. / order[i]) };
      if (mres < ires) mres=ires;
    }
  }
  return {acc.need_bulk, acc.need_deform, std::size_t(ceil(mres)),
          nsub_tidal_default, wdiv_tidal_default, acc.acc_radius};
}
  
namespace {
auto rho_from_dnu_switch(real_t dnu_switch, real_t gm1_center, 
                           const eos_barotr& eos) -> real_t
{
  const real_t gm1_switch{ 
    gm1_center + (1.0 + gm1_center) * std::expm1(-dnu_switch) 
  };
  assert(gm1_switch > 0);
  
  return eos.at_gm1(gm1_switch).rho();
}
}

auto tov_solver_fixstep::engine::get_deform(
            const eos_barotr& eos, 
            const spherical_star_info& prop,
            const std::vector<real_t>& dnu, 
            const std::vector<real_t>& rsqr, 
            const std::vector<real_t>& lambda, 
            const parameters& par) 
-> spherical_star_tidal
{
    if ((par.wdiv_tidal <= 0) || (par.wdiv_tidal >= 1)) {
      throw(std::runtime_error("get_deform: wdiv parameter"
                               "outside allowed range (0,1)"));
    }
    
    const std::size_t n1{ 
      std::max<std::size_t>(
                 5, std::size_t(dnu.size() * (1.-par.wdiv_tidal))) 
    };
    
    const std::size_t n2{ std::max<std::size_t>(5, dnu.size()-n1-1) };
    
    
    const real_t dnu_switch{ dnu[n1] };
    
    const real_t rho_switch{ 
      rho_from_dnu_switch(dnu_switch, prop.center_gm1, eos)
    };
    
    tidal_ode tode(eos, prop, dnu, rsqr, lambda, rho_switch);

    auto rtid{ integrate_ode_fixed(tode, par.nsub_tidal * n1) };

    real_t z_switch{ rtid[tidal_ode::YM2] };
    
    tidal_ode2 tode2( eos, prop, dnu, rsqr, lambda, 
                      rho_switch, z_switch);

    auto rtid2{ integrate_ode_fixed(tode2, par.nsub_tidal * n2) };

    return tode2.deformability(rtid2);
}
  

auto tov_solver_fixstep::engine::get_star_properties(
                      const eos_barotr& eos, 
                      const real_t rho_center, 
                      const parameters& par) 
-> spherical_star_properties  
{
  tov_ode ode{rho_center, eos, 1.001/real_t(par.nsamp_tov)};
  tov_ode::observer obs{ode};
  auto surf{ integrate_ode_fixed(ode, par.nsamp_tov, obs) };
  assert(obs.dnu.size()>0);

  auto prop{ ode.star(surf) };
  
  boost::optional<spherical_star_tidal> deform;
  if (par.find_tidal) {
    deform = get_deform(eos, prop, obs.dnu, 
                         obs.rsqr, obs.lambda, par);
  }

  boost::optional<spherical_star_bulk> bulk;
  if (par.find_bulk) {
    details::tov_profile prof{eos, prop, 
                obs.rsqr, obs.dnu, obs.lambda,
                obs.ebnd_by_r, obs.pvol_by_r};

    bulk = find_bulk_props(prof, par.bulk_acc); 
  }

  return {eos, prop, deform, bulk};  
}
  
  

auto tov_solver_fixstep::engine::get_star(
                      const eos_barotr& eos, 
                      const real_t rho_center,  
                      const parameters& par)
-> spherical_star
{
  tov_ode ode{rho_center, eos, 1.001/real_t(par.nsamp_tov)};
  tov_ode::observer obs{ode};
  auto surf{ integrate_ode_fixed(ode, par.nsamp_tov, obs) };
  assert(obs.dnu.size()>0);

  auto prop{ ode.star(surf) };
  
  boost::optional<spherical_star_tidal> deform;
  if (par.find_tidal) {
    deform = get_deform(eos, prop, obs.dnu, 
                         obs.rsqr, obs.lambda, par);
  }

  
  auto prof = std::make_shared<details::tov_profile>(eos, prop, 
                obs.rsqr, obs.dnu, obs.lambda,
                obs.ebnd_by_r, obs.pvol_by_r);

  boost::optional<spherical_star_bulk> bulk;
  if (par.find_bulk) {
    bulk = find_bulk_props(*prof, par.bulk_acc); 
  }

  return {prop, deform, bulk, prof};
}  
  
  
  

auto tov_solver_adaptive::engine::get_star_properties(
                      const eos_barotr& eos, 
                      const real_t rho_center, 
                      const parameters& par) 
-> spherical_star_properties
{
  tov_ode ode{rho_center, eos};
  tov_ode::observer obs{ode};
  auto surf{ 
    integrate_ode_adaptive(ode, par.acc_tov, par.nsamp_tov, obs) 
  };
  assert(obs.dnu.size()>0);

  auto prop{ ode.star(surf) };
  
  boost::optional<spherical_star_tidal> deform;
  if (par.find_tidal) {
    deform = get_deform(eos, prop, obs.dnu, 
                         obs.rsqr, obs.lambda, par);
  }

  boost::optional<spherical_star_bulk> bulk;
  if (par.find_bulk) {
    details::tov_profile prof{eos, prop, 
                obs.rsqr, obs.dnu, obs.lambda,
                obs.ebnd_by_r, obs.pvol_by_r};

    bulk = find_bulk_props(prof, par.bulk_acc); 
  }

  return {eos, prop, deform, bulk};  
}
  

auto tov_solver_adaptive::engine::get_star(
                      const eos_barotr& eos, 
                      const real_t rho_center,  
                      const parameters& par)
-> spherical_star
{
  tov_ode ode{rho_center, eos};
  tov_ode::observer obs{ode};
  auto surf{ 
    integrate_ode_adaptive(ode, par.acc_tov, par.nsamp_tov, obs) 
  };
  assert(obs.dnu.size()>0);

  auto prop{ ode.star(surf) };
  
  boost::optional<spherical_star_tidal> deform;
  if (par.find_tidal) {
    deform = get_deform(eos, prop, obs.dnu, 
                         obs.rsqr, obs.lambda, par);
  }

  
  auto prof = std::make_shared<details::tov_profile>(eos, prop, 
                obs.rsqr, obs.dnu, obs.lambda,
                obs.ebnd_by_r, obs.pvol_by_r);

  boost::optional<spherical_star_bulk> bulk;
  if (par.find_bulk) {
    bulk = find_bulk_props(*prof, par.bulk_acc); 
  }

  return {prop, deform, bulk, prof};
}    
  


auto tov_solver_adaptive::engine::get_deform(
            const eos_barotr& eos, 
            const spherical_star_info& prop,
            const std::vector<real_t>& dnu, 
            const std::vector<real_t>& rsqr, 
            const std::vector<real_t>& lambda, 
            const parameters& par)
-> spherical_star_tidal
{
    if ((par.wdiv_tidal <= 0) || (par.wdiv_tidal >= 1)) {
      throw(std::runtime_error("get_deform: wdiv parameter"
                               "outside allowed range (0,1)"));
    }
    
    const std::size_t n1{ 
      std::max<std::size_t>(
                 5, std::size_t(dnu.size() * (1.-par.wdiv_tidal))) 
    };
        
    const real_t dnu_switch{ dnu[n1] };
    
    const real_t rho_switch{ 
      rho_from_dnu_switch(dnu_switch, prop.center_gm1, eos)
    };

    
    tidal_ode tode(eos, prop, dnu, rsqr, lambda, rho_switch);

    auto rtid{ 
      integrate_ode_adaptive(tode, par.acc_tidal, par.nsamp_tidal) 
    };

    real_t z_switch{ rtid[tidal_ode::YM2]};
    
    tidal_ode2 tode2( eos, prop, dnu, rsqr, lambda, 
                      rho_switch, z_switch);

    auto rtid2{ 
      integrate_ode_adaptive(tode2, par.acc_tidal, par.nsamp_tidal) 
    };

    return tode2.deformability(rtid2);
}
  

auto get_tov_properties_adaptive(const eos_barotr eos, 
                   const real_t rho_center, const star_accuracy_spec acc) 
-> spherical_star_properties
{
  return tov_solver_adaptive::get_star_properties(eos, rho_center, acc);
}

auto get_tov_properties_adaptive(const eos_barotr eos, 
                   const real_t rho_center, 
                   const std::size_t nsamp_tov, 
                   const real_t acc_tov,
                   const std::size_t nsamp_tidal, 
                   const real_t acc_tidal,
                   const real_t wdiv_tidal,
                   const bool find_bulk, const bool find_tidal,
                   const real_t bulk_acc) 
-> spherical_star_properties
{
  return tov_solver_adaptive::engine::get_star_properties(
            eos, rho_center, {find_bulk, find_tidal, nsamp_tov,
            acc_tov, nsamp_tidal, acc_tidal, wdiv_tidal, bulk_acc});
}


auto get_tov_properties_fixstep(const eos_barotr eos, 
                   const real_t rho_center, const star_accuracy_spec acc) 
-> spherical_star_properties
{
  return tov_solver_fixstep::get_star_properties(eos, rho_center, acc);
}

auto get_tov_properties_fixstep(const eos_barotr eos, 
                   const real_t rho_center, 
                   const bool find_bulk, const bool find_tidal,
                   const std::size_t nsamp_tov, 
                   const std::size_t nsub_tidal,
                   const real_t wdiv_tidal,
                   const real_t bulk_acc)
-> spherical_star_properties
{
  return tov_solver_fixstep::engine::get_star_properties(
            eos, rho_center, { find_bulk, find_tidal, nsamp_tov, 
            nsub_tidal, wdiv_tidal, bulk_acc });
}

auto get_tov_properties(const eos_barotr eos, 
                   const real_t rho_center, 
                   const star_accuracy_spec acc) 
-> spherical_star_properties
{
  return get_tov_properties_adaptive(eos, rho_center, acc);
}

auto get_tov_star(const eos_barotr eos, 
                   const real_t rho_center, 
                   const star_accuracy_spec acc) 
-> spherical_star
{
  return tov_solver_adaptive::get_star(eos, rho_center, acc);
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



 
