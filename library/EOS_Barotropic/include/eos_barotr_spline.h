#ifndef EOS_BAROTR_SPLINE_H
#define EOS_BAROTR_SPLINE_H

#include "eos_barotropic.h"
#include <vector>
#include <functional>

namespace EOS_Toolkit {

/**\brief Create barotropic EOS based on splines

This EOS type internally uses regularly-spaced cardinal cubic splines.
The regular spacing implies that the evaluation cost is independent
from the table resolution.
Most quantities are interpolated with respect to logarithmic 
pseudo-enthalpy. In addition there is an interpolation spline 
that computes logarithmic pseudo-enthalpy from logarithmic density, 
which is used when the EOS is evaluated based on density.
The soundspeed is interpolated as function of logarithmic density.

The samples used for creating the EOS need to provide both 
pseudo-enthalpy and density, but may be arbitrarily spaced. 
They are resampled to the required spacing using
monotonic spline interpolation.

The range of the internal interpolation spline of the EOS has to be 
specified. It can span less than the range of the provided points. 
Below the interpolation range, the EOS will use a matching generalised 
polytropic EOS with the specified adiabatic index.
The density range that needs to be covered by the interpolation range
naturally depends on the application. It might therefore make sense to 
not use the full density range provided by the original EOS source.

One can specify whether the EOS is isentropic. If the temperature is 
zero, the EOS must be specified as isentropic or an exception is thrown.
Note that physically, the pseudo-enthalpy is already determined by 
pressure, density, and specific energy. For isentropic EOS, energy
density follows from density and pressure. The consistency of the 
provided samples is not checked and in the responsibility of the user. 
There are however similar functions which computes the pseudo-enthalpy 
itself using numerical integration, or both pseudo enthalpy and specific
energy (only possible for isentropic case).

The parameters are w.r.t EOS units, which have to be provided as a 
parameter and will be stored with the EOS for bookkeeping.


@param gm1 Samples for pseudo enthalpy \f$ g-1 \f$. Must be strictly 
           increasing.
@param rho Samples for mass density \f$ \rho \f$. Must be strictly 
           increasing and strictly positive.
@param eps Samples for specific energy \f$ \epsilon \f$. Must be 
           strictly increasing and obey \f$ \epsilon > -1 \f$
@param press Samples for pressure \f$ P \ge 0 \f$.
@param csnd Samples for soundspeed. Must obey 
           \f$ 0 \le c_s < 1 \f$
@param temp Samples for temperature or empty vector (equivalent 
            to all-zero). Must obey \f$ T\ge 0 \f$.
@param efrac Electron fraction or empty vector (in which case the EOS
             will not provide electron fraction).
@param isentropic Whether the EOS is isentropic, e.g. for degenerate 
                    matter.
@param rg_rho Range of \f$ \rho \f$ where EOS is computed from sample
                 points. Below, a polytrope is used. Range must be 
                 within the provided sample points.
@param n_poly Polytropic index used to extent EOS to zero density.
@param u   Unit system (w.r.t. SI) of the EOS
@param pts_per_mag Sample points per magnitude to be used internally
                   for EOS splines
                    
@return Generic interface employing tabulated barotropic EOS
**/
auto make_eos_barotr_spline(
  const std::vector<real_t>& gm1,
  const std::vector<real_t>& rho,
  const std::vector<real_t>& eps,
  const std::vector<real_t>& press,
  const std::vector<real_t>& csnd,
  const std::vector<real_t>& temp,
  const std::vector<real_t>& efrac, 
  bool isentropic,            
  interval<real_t> rg_rho,
  real_t n_poly,
  units u=units::geom_solar(),
  std::size_t pts_per_mag=200) 
-> eos_barotr;


/**\brief Create barotropic EOS based on splines

This EOS type internally uses regularly-spaced cardinal cubic splines.
The regular spacing implies that the evaluation cost is independent
from the table resolution.
Most quantities are interpolated with respect to logarithmic 
pseudo-enthalpy. In addition there is an interpolation spline 
that computes logarithmic pseudo-enthalpy from logarithmic density, 
which is used when the EOS is evaluated based on density.
The soundspeed is interpolated as function of logarithmic density.


The sample points provided for creating the EOS do not have to be
regularly spaced. They are resampled to the required spacing using
(arbitrary spaced) monotonic spline interpolation.

The samples are provided with respect to density. The pseudo-enthalpy
does not have to be provided, but is computed using numerical 
integration via the trapezoidal rule.
If sample points of density are further apart than the specified 
minimum points per magnitude, additional sample points are inserted
via cubic monotonic spline interpolation before the numerical 
integration. This ensures that the final EOS is a consistent 
representation of a cubic monotonic spline interpolation of the 
sample points with respect to density.

Note there is a similar function where the pseudo-enthalpy is provided
as samples as well. This should be used if it is know accurately, e.g
for EOS based on analytic expressions. Another function specially for 
isentropic EOS recomputes the specific energy as well as the pseudo 
enthalpy.

The range of the internal interpolation spline of the final EOS has to 
be specified. It can span less than the range of the provided points. 
Below the interpolation range, the EOS will use a matching generalised 
polytropic EOS with the specified adiabatic index.
The density range that needs to be covered by the interpolation range
naturally depends on the application. It might therefore make sense to 
not use the full density range provided by the original EOS source.

One can specify whether the EOS is isentropic. If the temperature is 
zero, the EOS must be specified as isentropic or an exception is thrown.
Note that physically, energy density follows from density and pressure 
for the case of isentropic EOS. The consistency of the provided samples 
is not checked and in the responsibility of the user. 

The parameters are w.r.t EOS units, which have to be provided as a 
parameter and will be stored with the EOS for bookkeeping.


@param rho Samples for mass density \f$ \rho \f$. Must be strictly 
           increasing and strictly positive.
@param eps Samples for specific energy \f$ \epsilon \f$. Must be 
           strictly increasing and obey \f$ \epsilon > -1 \f$
@param press Samples for pressure \f$ P \ge 0 \f$.
@param csnd Samples for soundspeed. Must obey 
           \f$ 0 \le c_s < 1 \f$
@param temp Samples for temperature or empty vector (equivalent 
            to all-zero). Must obey \f$ T\ge 0 \f$.
@param efrac Electron fraction or empty vector (in which case the EOS
             will not provide electron fraction).
@param isentropic Whether the EOS is isentropic, e.g. for degenerate 
                    matter.
@param rg_rho Range of \f$ \rho \f$ where EOS is computed from sample
                 points. Below, a polytrope is used. Range must be 
                 within the provided sample points.
@param n_poly Polytropic index used to extent EOS to zero density.
@param u   Unit system (w.r.t. SI) of the EOS
@param pts_per_mag Sample points per magnitude to be used internally
                   for EOS splines
                    
@return Generic interface employing tabulated barotropic EOS
**/
auto make_eos_barotr_spline(
  const std::vector<real_t>& rho,
  const std::vector<real_t>& eps,
  const std::vector<real_t>& press,
  const std::vector<real_t>& csnd,
  const std::vector<real_t>& temp,
  const std::vector<real_t>& efrac, 
  bool isentropic,            
  interval<real_t> rg_rho,
  real_t n_poly,
  units u=units::geom_solar(),
  std::size_t pts_per_mag=200)
->eos_barotr;


/**\brief Create barotropic EOS based on splines

This EOS type internally uses regularly-spaced cardinal cubic splines.
The regular spacing implies that the evaluation cost is independent
from the table resolution.
Most quantities are interpolated with respect to logarithmic 
pseudo-enthalpy. In addition there is an interpolation spline 
that computes logarithmic pseudo-enthalpy from logarithmic density, 
which is used when the EOS is evaluated based on density.
The soundspeed is interpolated as function of logarithmic density.

The sample points provided for creating the EOS do not have to be
regularly spaced. They are resampled to the required spacing using
(arbitrary spaced) monotonic spline interpolation.

The samples are provided with respect to density. The pseudo-enthalpy
does not have to be provided, but is computed using numerical 
integration via the trapezoidal rule. It is assumed that the EOS
is isentropic. The specific energy is recomputed from the pressure
and density samples under this assumption.
If sample points of density are further apart than the specified 
minimum points per magnitude, additional sample points are inserted
via cubic monotonic spline interpolation before the numerical 
integrations of pseudo-enthalpy and specific energy. This ensures 
that the final EOS is a consistent representation of a cubic monotonic 
spline interpolation of the sample points with respect to density.

Note there are similar functions where specific energy or specific 
energy and pseudo-enthalpy are provided as samples as well. This should 
be used if they are know accurately, e.g, for EOS based on analytic 
expressions.

The range of the internal interpolation spline of the final EOS has to 
be specified. It can span less than the range of the provided points. 
Below the interpolation range, the EOS will use a matching generalised 
polytropic EOS with the specified adiabatic index.
The density range that needs to be covered by the interpolation range
naturally depends on the application. It might therefore make sense to 
not use the full density range provided by the original EOS source.

If a temperature is given, its consistency with the assumption of 
isentropy cannot be checked and is in the responsability of the user.

The parameters are w.r.t EOS units, which have to be provided as a 
parameter and will be stored with the EOS for bookkeeping.

@param rho Samples for mass density \f$ \rho \f$. Must be strictly 
           increasing and strictly positive.
@param press Samples for pressure \f$ P \ge 0 \f$.
@param csnd Samples for soundspeed. Must obey 
           \f$ 0 \le c_s < 1 \f$
@param temp Samples for temperature or empty vector (equivalent 
            to all-zero). Must obey \f$ T\ge 0 \f$.
@param efrac Electron fraction or empty vector (in which case the EOS
             will not provide electron fraction).
@param rg_rho Range of \f$ \rho \f$ where EOS is computed from sample
                 points. Below, a polytrope is used. Range must be 
                 within the provided sample points.
@param n_poly Polytropic index used to extent EOS to zero density.
@param eps0    Specific energy at zero density
@param u   Unit system (w.r.t. SI) of the EOS
@param pts_per_mag Sample points per magnitude to be used internally
                   for EOS splines
                    
@return Generic interface employing tabulated barotropic EOS
**/
auto make_eos_barotr_spline(
  const std::vector<real_t>& rho,
  const std::vector<real_t>& press,
  const std::vector<real_t>& csnd,
  const std::vector<real_t>& temp,
  const std::vector<real_t>& efrac, 
  interval<real_t> rg_rho,
  real_t n_poly,
  real_t eps0,
  units u=units::geom_solar(),
  std::size_t pts_per_mag=200)
-> eos_barotr;


auto make_eos_barotr_spline(
  std::function<real_t(real_t)> gm1_rho, 
  std::function<real_t(real_t)> rho_gm1, 
  std::function<real_t(real_t)> eps_gm1, 
  std::function<real_t(real_t)> press_gm1, 
  std::function<real_t(real_t)> csnd_rho, 
  std::function<real_t(real_t)> temp_gm1, 
  std::function<real_t(real_t)> efrac_gm1, 
  bool isentropic, 
  interval<real_t> rg_rho, 
  real_t n_poly,
  units u, 
  std::size_t pts_per_mag=200)
->eos_barotr;


/**\brief Create barotropic EOS based on splines

This allows to represent existing EOS objects of any type 
as a spline EOS by sampling them
it.

@param eos The EOS to be sampled
@param rg_rho The density range that should be represented by
              interpolating splines
@param n_poly Polytropic index used to extent EOS to zero density.
@param pts_per_mag Sample points per magnitude to be used internally
                   for EOS splines
                    
@return Generic interface employing tabulated barotropic EOS
*/
auto make_eos_barotr_spline(
  const eos_barotr& eos, 
  interval<real_t> rg_rho, 
  real_t n_poly,
  std::size_t pts_per_mag=200)
-> eos_barotr;

}//namespace EOS_Toolkit


#endif

