#ifndef EOS_BAROTR_SPLINE_IMPL_H
#define EOS_BAROTR_SPLINE_IMPL_H

#include "eos_barotropic_impl.h"
#include "eos_barotr_gpoly_impl.h"
#include "interpol_logspl.h"
#include <vector>
#include <boost/optional.hpp>

namespace EOS_Toolkit {
namespace implementations {


///Barotropic EOS based on cardinal spline interpolation. 
/**
This uses  logarithmically spaced cardinal spline interpolation.
For densities below the smallest tabulated value, 
a matching polytropic EOS (see eos_genpoly) is used. 
For notation, see eos_cold. 
*/
class eos_barotr_spline : public eos_barotr_impl {
   
  public:
  
  using lgspl_t = detail::interpol_logspl_impl;
  using lglgspl_t = detail::interpol_llogspl_impl;
  using opt_t = boost::optional<lgspl_t>;
  

  ///Constructor
  eos_barotr_spline(
    lglgspl_t gm1_,   ///< \f$ g - 1 \f$ from \f$ \rho \f$
    lglgspl_t rho_,   ///< \f$ \rho \f$ from \f$ g - 1 \f$
    lgspl_t eps_,   ///< \f$ \epsilon \f$ from \f$ g - 1 \f$
    lglgspl_t p_,   ///< \f$ P \f$ from \f$ g - 1 \f$
    lgspl_t hm1_,   ///< \f$ \epsilon + P/\rho \f$ from \f$ g - 1 \f$
    lgspl_t csnd_,   ///< \f$ c_s^2 \f$ from \f$ g - 1 \f$
    opt_t temp_,  ///< \f$ T \f$ from \f$ g - 1 \f$ (optional) 
    opt_t efrac_, ///< \f$ Y_e \f$ from \f$ g - 1 \f$ (optional)
    bool isentropic_,  ///< Whether EOS is isentropic
    const eos_barotr_gpoly& poly_ ///< Polytropic EOS for low densities
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
  bool has_efrac() const final  {return bool(efrac_gm1);}  
  
  
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

  ///Compute adiabatic soundspeed \f$ c_s \f$
  /**Assumes input is in the valid range, no checks are performed.*/
  real_t csnd_from_rho_gm1(
    real_t rho,     ///<Rest mass density  \f$ \rho \f$
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

  void save(datasink s) const final;
  auto descr_str() const -> std::string final;

  const static bool file_handler_registered;
  static const std::string datastore_id;

  private:

  static auto get_rggm1(const lgspl_t&, const lglgspl_t&, 
                        const lgspl_t&, const lglgspl_t&, 
                        const opt_t&, const opt_t&)  
  -> range;

  lglgspl_t gm1_rho;
  lgspl_t eps_gm1; 
  lglgspl_t p_gm1;
  lgspl_t hm1_gm1;
  lglgspl_t rho_gm1;
  lgspl_t csnd_rho;
  opt_t temp_gm1;
  opt_t efrac_gm1;

  const eos_barotr_gpoly poly;

  const range rggm1;
  const range rgrho;
  const real_t gm1_low;
  const real_t rho_low;
  const real_t min_h{1.};
  real_t efrac0{0.};
  real_t temp0{0.};
  
  bool zerotemp{true};    ///< If EOS is zero temperature
  const bool isentropic;  ///< If EOS is isentropic
};

}//namespace implementations 

}//namespace EOS_Toolkit

#endif

