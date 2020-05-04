#ifndef EOS_BAROTR_TABLE_IMPL_H
#define EOS_BAROTR_TABLE_IMPL_H

#include "eos_barotropic_impl.h"
#include "eos_barotr_gpoly_impl.h"
#include "interpol.h"
#include <vector>

namespace EOS_Toolkit {
namespace implementations {


///Tabulated barotropic EOS. 
/**
This uses lookup tables with logarithmic interpolation (lookup_table).
The interpolation is not thermodynamically consistent, so
tables with a large number of points (>1000) should be used.
For densities below the smallest tabulated value, 
a matching polytropic EOS (see eos_genpoly) is used. 
For notation, see eos_cold. 
*/
class eos_barotr_table : public eos_barotr_impl {
  bool zerotemp{true};    ///< If EOS is zero temperature
  const bool isentropic;  ///< If EOS is isentropic
  const bool hasefrac;    ///< If EOS has electron fraction information

  const range rgrho;
  const range rggm1;

  lookup_table_magx gm1_rho, eps_gm1, hm1_gm1;
  lookup_table_magx pbr_gm1, rho_gm1, cs2_gm1;
  lookup_table_magx temp_gm1{};
  lookup_table_magx efrac_gm1{};

  real_t min_h;
  real_t efrac0{0.};
  real_t temp0{0.};
  
  eos_barotr_gpoly poly;
   
  public:

  using func_t  = std::function<real_t(real_t)>;
  

  ///Constructor
  eos_barotr_table(
    range rg_rho_, ///< \f$ \rho \f$ range to tabulate
    range rg_gm1_, ///< \f$ g-1 \f$ range to tabulate
    std::size_t nsamples_, ///< Number of samples for lookup table
    int magnitudes_,       ///< Dynamic range for lookup table
    func_t gm1_,   ///< \f$ g - 1 \f$ from \f$ \rho \f$
    func_t rho_,   ///< \f$ \rho \f$ from \f$ g - 1 \f$
    func_t eps_,   ///< \f$ \epsilon \f$ from \f$ g - 1 \f$
    func_t pbr_,   ///< \f$ P/\rho \f$ from \f$ g - 1 \f$
    func_t cs2_,   ///< \f$ c_s^2 \f$ from \f$ g - 1 \f$
    func_t temp_,  ///< \f$ T \f$ from \f$ g - 1 \f$ (or nullptr) 
    func_t efrac_, ///< \f$ Y_e \f$ from \f$ g - 1 \f$ (or nullptr)
    bool isentropic_,  ///< Whether EOS is isentropic
    const eos_barotr_gpoly& poly_  ///< Polytropic EOS for low densities
  );

  
  ///Returns range of validity for density
  const range& range_rho() const final {return rgrho;}
  
  ///Returns range of validity for \f$ g-1 \f$ 
  const range& range_gm1() const final {return rggm1;}

  ///Returns range of validity for \f$ g-1 \f$ 
  real_t minimal_h() const final {return min_h;}
  
  
  ///Whether EOS is isentropic
  bool is_isentropic() const final  {return isentropic;}
  
  ///Whether EOS is for zero temperature
  bool is_zero_temp() const final  {return zerotemp;}
  
  ///Whether EOS can compute temperature
  bool has_temp() const final  {return true;}
  
  ///Whether EOS can compute electron fraction
  bool has_efrac() const final  {return hasefrac;}  
  
  
  ///Compute \f$ g-1 \f$ 
  /**Assumes input is in the valid range, no checks are performed.*/
  real_t gm1_from_rho(
    real_t rho      ///<Rest mass density  \f$ \rho \f$
  ) const final;

  ///Compute Rest mass density \f$ \rho \f$
  /**Assumes input is in the valid range, no checks are performed.*/
  real_t rho(
    real_t gm1      ///< \f$ g-1 \f$
  ) const final;


  ///Compute Specific internal energy \f$\epsilon \f$
  /**Assumes input is in the valid range, no checks are performed.*/
  real_t eps(
    real_t gm1      ///< \f$ g-1 \f$
  ) const final; 


  ///Compute Pressure \f$ P \f$ 
  /**Assumes input is in the valid range, no checks are performed.*/
  real_t press(
    real_t gm1      ///< \f$ g-1 \f$
  ) const final;


  ///Compute specific enthalpy (excluding restmass) \f$ h-1  \f$
  /**Assumes input is in the valid range, no checks are performed.*/
  real_t hm1(
    real_t gm1      ///< \f$ g-1 \f$
  ) const final;

  ///Compute adiabatic soundspeed \f$ c_s \f$
  /**Assumes input is in the valid range, no checks are performed.*/
  real_t csnd(
    real_t gm1      ///< \f$ g-1 \f$
  ) const final;

  ///Returns temperature \f$ T \f$
  /**Assumes input is in the valid range, no checks are performed.*/
  real_t temp(
    real_t gm1      ///< \f$ g-1 \f$
  ) const final;

  ///Compute electron fraction \f$ Y_e \f$ 
  /**Assumes input is in the valid range, no checks are performed.*/
  real_t ye(
    real_t gm1      ///< \f$ g-1 \f$
  ) const final;

};

}//namespace implementations 

}//namespace EOS_Toolkit

#endif

