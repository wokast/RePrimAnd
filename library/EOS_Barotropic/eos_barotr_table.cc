#include "eos_barotr_table.h"
#include "eos_barotr_table_impl.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>

using namespace std;
using namespace EOS_Toolkit;
using namespace EOS_Toolkit::implementations;
//~ using namespace EOS_Toolkit::detail;


/**

*/
eos_barotr_table::eos_barotr_table(
        range rg_rho_, range rg_gm1_, 
        std::size_t nsamples_, int magnitudes_,
        func_t gm1_, func_t rho_, func_t eps_,  func_t pbr_,   
        func_t cs2_, func_t temp_, func_t efrac_, 
        bool isentropic_, const eos_barotr_gpoly& poly_)
: isentropic{isentropic_}, hasefrac{efrac_},
  rgrho{0, rg_rho_.max()},
  rggm1{0, rg_gm1_.max()},
  gm1_rho{std::move(gm1_), rg_rho_ , nsamples_, magnitudes_}, 
  eps_gm1{std::move(eps_), rg_gm1_, nsamples_, magnitudes_},
  pbr_gm1{std::move(pbr_), rg_gm1_, nsamples_, magnitudes_}, 
  rho_gm1{std::move(rho_), rg_gm1_, nsamples_, magnitudes_}, 
  cs2_gm1{std::move(cs2_), rg_gm1_, nsamples_, magnitudes_}, 
  poly{poly_}
{
  if ((rho_gm1.range_y().min() < 0.0) 
        || (gm1_rho.range_x().min() < 0.0)) {
    throw runtime_error("eos_barotr_table: negative mass density");
  }
  if (cs2_gm1.range_y().max() >= 1.0) {
    throw runtime_error("eos_barotr_table: superluminal sound speed");
  }
  if (cs2_gm1.range_y().min() <= 0.0) {
    throw runtime_error("eos_barotr_table: imaginary sound speed");
  }
  if (pbr_gm1.range_y().min() <= 0.0) {
    throw runtime_error("eos_barotr_table: pressure <= 0");
  }
  if (gm1_rho.range_y().min() < 0.0) {
    throw runtime_error("eos_barotr_table: encountered g < 1");
  }
  
  if (temp_) {
    temp_gm1 = {std::move(temp_), rg_gm1_, nsamples_, magnitudes_};
    temp0    = temp_gm1(rg_gm1_.min());
    if (temp_gm1.range_y().min() < 0.0) {
      throw runtime_error("eos_barotr_table: encountered negative "
                          "temperature");
    }
    zerotemp = temp_gm1.range_y().max() == 0;
  }
  
  if (zerotemp && !isentropic) {
    throw runtime_error("eos_barotr_table: zero-temperature EOS must "
                        "be isentropic");
  }
  
  if (hasefrac) {
    efrac_gm1 = {std::move(efrac_), rg_gm1_, nsamples_, magnitudes_};
    efrac0    = efrac_gm1(rg_gm1_.min());
  }
  
  auto h1g1 = [this] (real_t g1) {
    return eps_gm1(g1) + pbr_gm1(g1);
  };
  hm1_gm1   = {h1g1, rg_gm1_, nsamples_, magnitudes_};
  min_h     = 1.0 + std::min(poly.hm1(0.), hm1_gm1.range_y().min());
}

real_t eos_barotr_table::gm1_from_rho(real_t rho) const
{
  return (rho > gm1_rho.range_x().min()) 
            ? gm1_rho(rho) : poly.gm1_from_rho(rho); 
}


real_t eos_barotr_table::eps(real_t gm1) const
{
  return (gm1 > eps_gm1.range_x().min()) 
            ? eps_gm1(gm1) : poly.eps(gm1);
}


real_t eos_barotr_table::press(real_t gm1) const
{
  return (gm1 > pbr_gm1.range_x().min()) 
          ? (pbr_gm1(gm1)*rho_gm1(gm1)) : poly.press(gm1);
}

real_t eos_barotr_table::rho(real_t gm1) const
{
  return (gm1 > rho_gm1.range_x().min()) 
            ? rho_gm1(gm1) : poly.rho(gm1);
}

real_t eos_barotr_table::hm1(real_t gm1) const
{
  return (gm1 > hm1_gm1.range_x().min()) 
            ? hm1_gm1(gm1)  : poly.hm1(gm1);
}

real_t eos_barotr_table::csnd(real_t gm1) const
{
  return (gm1 > cs2_gm1.range_x().min()) 
            ? sqrt(cs2_gm1(gm1)) : poly.csnd(gm1);
}

real_t eos_barotr_table::temp(real_t gm1) const
{
  if (zerotemp) return 0;
  return (gm1 < temp_gm1.range_x()) ? temp0 : temp_gm1(gm1);
}

real_t eos_barotr_table::ye(real_t gm1) const
{
  if (!has_efrac())
    throw runtime_error("eos_barotr_table: electron fraction "
                        "not available.");
  return (gm1 > efrac_gm1.range_x().min()) ? efrac_gm1(gm1) : efrac0;
}



eos_barotr EOS_Toolkit::make_eos_barotr_table(
  const std::vector<real_t>& gm1, const std::vector<real_t>& rho,
  const std::vector<real_t>& eps, const std::vector<real_t>& pbr, 
  const std::vector<real_t>& cs2, const std::vector<real_t>& temp, 
  const std::vector<real_t>& efrac, bool isentropic, real_t n_poly)
{
  const size_t tsize = rho.size();
  if (tsize<5) {
    throw std::invalid_argument("make_eos_barotr_table: want at least "
                                "5 sample points");
  }
  if ((tsize != eps.size()) || (tsize != pbr.size()) ||  
      (tsize != cs2.size()) || (tsize != gm1.size()) ||
      ((!temp.empty()) && (temp.size() != tsize)) ||
      ((!efrac.empty()) && (efrac.size() != tsize))) 
  {
    throw std::invalid_argument("make_eos_barotr_table: number of "
                      "samples for different quantities don't match");
  }
  if (rho[0]<=0) {
    throw std::invalid_argument("make_eos_barotr_table: sampling "
                                "densities must be strictly positive");
  }
  
  auto poly = eos_barotr_gpoly::from_boundary(rho[0], eps[0], 
                                      pbr[0]*rho[0], n_poly, rho[1]);
                                       
  const real_t gsc = (poly.gm1_from_rho(rho[0]) - gm1[0]) 
                      / (1.0 + gm1[0]);
  
  std::vector<real_t> ngm1;
  for (auto g1 : gm1) ngm1.push_back(g1 + gsc * (1.0 + g1));
  
  assert(ngm1[0]>0);
  
  cspline_mono sgm1{rho, ngm1};
  cspline_mono srho{ngm1, rho};
  cspline_mono seps{ngm1, eps};
  cspline_mono spbr{ngm1, pbr};
  cspline_mono scs2{ngm1, cs2};
  
  eos_barotr_table::func_t stemp{nullptr};
  if (!temp.empty()) stemp = cspline_mono{ngm1, temp};
  
  eos_barotr_table::func_t sefrac{nullptr};
  if (!efrac.empty()) sefrac = cspline_mono{ngm1, efrac};
  
  const std::size_t nsample = 10*rho.size(); //TODO: better criterion
  
  const int i50       = int(ceil(rho.size()/2.));
  const real_t rho50 = rho.at(i50);
  if ((rho50 <= 0) || (rho.back() <= rho50)) {
    throw std::invalid_argument("make_eos_barotr_table: sample "
                               "points for density invalid");
  }
  const int mags = max(1, int(2*log10(rho.back() / rho50)));
  
  
  return eos_barotr{ std::make_shared<eos_barotr_table>(
    sgm1.range_x(), srho.range_x(), nsample, mags,
    sgm1, srho, seps, spbr, scs2, stemp, sefrac, isentropic, poly) 
  };
}

