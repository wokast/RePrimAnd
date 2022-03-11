#ifndef SPHERICAL_STARS_H
#define SPHERICAL_STARS_H
#include <memory>
#include <boost/optional.hpp>
#include "config.h"
#include "eos_barotropic.h"

namespace EOS_Toolkit {

/// Spherical star tidal deformability measures
struct spherical_star_tidal {
  real_t k2;     ///< Love number \f$ k_2 \f$ (dimensionless)
  real_t lambda; ///< Dimensionless tidal deformability \f$ \Lambda \f$
};

/// Spherical star bulk properties
struct spherical_star_bulk {
  real_t circ_radius;    ///< Circumferential radius 
  real_t rho;            ///< Baryonic mass density
  real_t proper_volume;  ///< Enclosed proper volume
  real_t bary_mass;      ///< Enclosed baryonic mass
};

/// Spherical star assorted measures
struct spherical_star_info {
  real_t center_rho;       ///< Central baryonic mass density
  real_t center_gm1;       ///< Central pseudo enthalpy \f$ g - 1 \f$
  real_t center_nu;        ///< Central metric potential \f$ \nu \f$
  real_t grav_mass;        ///< Gravitational (AMD) mass
  real_t binding_energy;   ///< Binding energy
  real_t circ_radius;      ///< Circumferential surface radius
  real_t proper_volume;    ///< Proper volume
  real_t moment_inertia;   ///< Moment of inertia
};

/**\brief Class to describe a spherical neutron star

This collects scalar measures of a spherical neutron star
solution, but not the radial profile. It also provides the EOS.
All quantities are in the same geometric units used by the EOS.
**/
class spherical_star_properties {
  public: 
  
  using deform_t = boost::optional<spherical_star_tidal>;
  using bulk_t = boost::optional<spherical_star_bulk>;
    
  ///Copy constructor.
  spherical_star_properties(
                     const spherical_star_properties &) = default;
  
  ///Move constructor
  spherical_star_properties(spherical_star_properties &&) = default;
  
  ///Assignment operator.
  spherical_star_properties& operator=(
                     const spherical_star_properties&) = default;
  
  ///Move assignment operator
  spherical_star_properties& operator=(
                     spherical_star_properties&&) = default;

  ///The constructor 
  spherical_star_properties(eos_barotr eos_, spherical_star_info info_,
                            deform_t deform_, bulk_t bulk_);
                            
  /// Central matter state
  auto center_state() const -> eos_barotr::state; 
  
  /// Central baryonic mass density
  auto center_rho() const -> real_t;
  
  /// Central pseudo enthalpy \f$ g - 1 \f$
  auto center_gm1() const -> real_t;
  
  /// Central specific internal energy \f$ \epsilon \f$
  auto center_eps() const -> real_t;
  
  /// Central pressure
  auto center_press() const -> real_t;
  
  /// Central sound speed
  auto center_csnd() const -> real_t;
  
  /// Central electron fraction
  /**
  This is only available if supported by the EOS, otherwise
  an exception is thrown.
  **/
  auto center_ye() const -> real_t;
  
  /// Gravitational mass \f$ M_g \f$
  auto grav_mass() const -> real_t;
  
  /// Baryonic mass \f$ M_b \f$
  auto bary_mass() const -> real_t;
  
  /// Binding energy \f$ M_b - M_g \f$
  auto binding_energy() const -> real_t;
  
  /// Circumferential surface radius
  auto circ_radius() const -> real_t;
  
  /// Proper volume
  auto proper_volume() const -> real_t;
  
  /// Moment of inertia
  auto moment_inertia() const -> real_t;

  /// Whether the tidal deformability is known
  auto has_deform() const -> bool;
  
  /// Whether the bulk properties are known
  auto has_bulk() const -> bool;
  
  /// Tidal deformability measures
  /**
  Throws exception if not available
  **/
  auto deformability() const -> const spherical_star_tidal&;
  
  /// Bulk properties
  /**
  Throws exception if not available
  **/  
  auto bulk() const -> const spherical_star_bulk&;

  /// The EOS
  auto eos() const -> const eos_barotr&;
  
  private: 
  
  const eos_barotr _eos;
  spherical_star_info _info;
  deform_t _deform;
  bulk_t  _bulk;
};

class spherical_star_profile {
  const eos_barotr _eos;
  real_t _surf_radius;
  
  public:
  explicit spherical_star_profile(eos_barotr eos_, real_t surf_radius_);
  virtual ~spherical_star_profile() = default;
  
  auto eos() const -> const eos_barotr&;
  
  auto surf_circ_radius() const -> real_t;

  virtual real_t center_gm1() const=0;
  virtual real_t nu_from_rc(real_t rc) const=0;
  virtual real_t lambda_from_rc(real_t rc) const=0;
  virtual real_t gm1_from_rc(real_t rc) const=0;  
  virtual real_t mbary_from_rc(real_t rc) const=0;  
  virtual real_t pvol_from_rc(real_t rc) const=0;  
  
  auto state_from_rc(real_t rc) const -> eos_barotr::state;

  //~ auto g_rcrc_from_rc(real_t rc) const -> real_t;  
  //~ auto g_phph_from_rc_th(real_t rc, real_t th) const -> real_t;  
  //~ auto g_thth_from_rc(real_t rc) const -> real_t;  
  //~ auto lapse_from_rc(real_t rc) const -> real_t;  
};

/**\brief Class representing a spherical neutron star model

This provides everything spherical_star_properties does,
and in addition the radial profiles for metric and matter 
quantities. The profiles are provided as functions of circumferential
radius that can be evaluated both inside and outside the star.
All quantities are in the same geometric units used by the EOS.
**/
class spherical_star : public spherical_star_properties  {
  using pprof_t = std::shared_ptr<const spherical_star_profile>;
  
  pprof_t pprof; 
  auto profile() const -> const spherical_star_profile&;

  public:
  
  ///Constructor.
  spherical_star(spherical_star_info info_,
                 deform_t deform_, bulk_t bulk_,
                 pprof_t pprof_);

 
  
  ///Copy constructor.
  spherical_star(const spherical_star &)  = default;
  
  ///Move constructor
  spherical_star(spherical_star &&)       = default;
  
  ///Assignment operator.
  spherical_star& operator=(const spherical_star&)  = default;
  
  ///Move assignment operator
  spherical_star& operator=(spherical_star&&)       = default;

  /// Metric potential \f$ \nu \f$ at given circumferential radius
  auto nu_from_rc(real_t rc) const -> real_t;
  
  /// Metric potential \f$ \lambda \f$ at given circumferential radius
  auto lambda_from_rc(real_t rc) const -> real_t;
  
  /// Baryonic mass within given circumferential radius
  auto mbary_from_rc(real_t rc) const -> real_t;
  
  /// Proper volume within given circumferential radius
  auto pvol_from_rc(real_t rc) const -> real_t;
  
  /// Matter state at circumferential radius
  auto state_from_rc(real_t rc) const -> eos_barotr::state;

  /// Pseudo enthalpy \f$ g - 1 \f$ at circumferential radius
  auto gm1_from_rc(real_t rc) const -> real_t;
  
  /// Baryonic mass density \f$ \rho \f$ at circumferential radius
  auto rho_from_rc(real_t rc) const -> real_t;  
  
  /// Pressure  \f$ P \f$ at circumferential radius
  auto press_from_rc(real_t rc) const -> real_t;  
  
  /// Specific internal energy \f$ \epsilon \f$ at circumferential radius
  auto eps_from_rc(real_t rc) const -> real_t;  
  
  /// Sound speed \f$ c_s \f$ at circumferential radius
  auto csnd_from_rc(real_t rc) const -> real_t;  
  
  /// Electron fraction \f$ Y_e \f$ at circumferential radius
  /**
  This is only available if supported by the EOS, otherwise
  an exception is thrown.
  **/
  auto ye_from_rc(real_t rc) const -> real_t;  
  
  /// Temperature \f$ T [MeV] \f$ at circumferential radius
  /**
  This is only available if supported by the EOS, otherwise
  an exception is thrown.
  **/
  auto temp_from_rc(real_t rc) const -> real_t;  
  
  
};


///Class for specifiying heuristic accuracy of TOV solution
struct tov_acc_simple {  
  const real_t tov;    ///< Accuracy for adaptive solution of TOV ODE
  const real_t deform; ///< Accuracy for adaptive solution of tidal ODE
  
  /// Minimum number of steps used for solving ODEs
  const std::size_t minsteps; 
  
  ///Constructor
  tov_acc_simple(real_t acc_tov_=1e-8, 
                 real_t acc_deform_=1e-6,
                 std::size_t minsteps_=500)
  : tov{acc_tov_}, deform{acc_deform_}, minsteps{minsteps_} {}
};

///Class for specifying absolute tolerances for TOV solution
struct tov_acc_precise {  
  const real_t mass;     ///< Tolerance for masses
  const real_t radius;   ///< Tolerance for radius
  const real_t minertia; ///< Tolerance for moment of inertia
  const real_t deform;   ///< Tolerance for  tidal deformability
  
  /// Minimum number of steps used for solving ODEs
  const std::size_t minsteps; 
  /// Minimum adaptive tolerance to try before giving up
  const real_t acc_min; 
  
  ///Constructor
  tov_acc_precise(real_t acc_mass_=1e-8, 
                  real_t acc_radius_=1e-8, 
                  real_t acc_minertia_=1e-8, 
                  real_t acc_deform_=1e-6,
                  std::size_t minsteps_=500,
                  real_t acc_min_=1e-14)
  : mass{acc_mass_}, radius{acc_radius_}, minertia{acc_minertia_},
    deform{acc_deform_}, minsteps{minsteps_}, acc_min{acc_min_} {}
};


/**\brief Compute spherical neutron star model. 

@param eos The (barotropic) EOS of the NS. 
@param rho_center The central baryonic mass density. Units are the 
same as used by the EOS.
@param acc Tolerances for adaptive ODE solver
@param find_bulk Option to also compute the "bulk" properties 
@param find_tidal Whether to compute the tidal deformability

@return Stellar model 

This returns a class providing all global NS properties as well as
the stellar profile. 
**/ 
auto make_tov_star(const eos_barotr eos, const real_t rho_center, 
                   const tov_acc_simple acc, 
                   const bool find_bulk=false, 
                   const bool find_tidal=true) 
-> spherical_star;


/**\brief Compute properties of spherical neutron star. 

@param eos The (barotropic) EOS of the NS. 
@param rho_center The central baryonic mass density. Units are the 
same as used by the EOS.
@param acc Tolerances for adaptive ODE solver
@param find_bulk Option to also compute the "bulk" properties 
@param find_tidal Whether to compute the tidal deformability

@return Stellar model properties

This returns a class providing all global NS properties but not
the stellar profile. 
**/ 
auto get_tov_star_properties(const eos_barotr eos, 
                   const real_t rho_center, const tov_acc_simple acc, 
                   const bool find_bulk=false, 
                   const bool find_tidal=true) 
-> spherical_star_properties;


/**\brief Compute spherical neutron star model. 

@param eos The (barotropic) EOS of the NS. 
@param rho_center The central baryonic mass density. Units are the 
same as used by the EOS.
@param acc Absolute error bounds for individual stellar properties
@param find_bulk Option to also compute the "bulk" properties 
@param find_tidal Whether to compute the tidal deformability

@return Stellar model 

This returns a class providing all global NS properties as well as
the stellar profile. 
**/ 
auto make_tov_star(const eos_barotr eos, const real_t rho_center, 
              const tov_acc_precise acc, const bool find_bulk=false, 
              const bool find_tidal=true)
-> spherical_star;


/**\brief Compute properties of spherical neutron star. 

@param eos The (barotropic) EOS of the NS 
@param rho_center The central baryonic mass density. Units are the 
same as used by the EOS
@param acc Absolute error bounds for individual stellar properties
@param find_bulk Option to also compute the "bulk" properties 
@param find_tidal Whether to compute the tidal deformability

@return Stellar model properties

This returns a class providing all global NS properties but not
the stellar profile. 
**/ 
auto get_tov_star_properties(const eos_barotr eos, 
                   const real_t rho_center, 
                   const tov_acc_precise acc, 
                   const bool find_bulk=false, 
                   const bool find_tidal=true)
-> spherical_star_properties;


/**\brief Find maximum mass TOV model
@param eos The (barotropic) EOS of the NS
@param rhobr0 Lower bound of central density for search
@param rhobr1 Upper bound of central density for search
@param bits Accuracy of maximum search (number significant bits)
@param acc Accuracy for TOV solving
@param max_steps Maximum steps for maximum search

@return Central baryonic mass density of maximum mass model

This finds the model with maximum gravitational mass within
a specified range of the central baryonic mass density. The
provided range is first limited to the validity range of the EOS.
Since the EOS framework guarantees causality, the returned value
might be the one for the maximum causal central density, instead
of a true maximum (i.e. zero derivative). Similarly, if the 
provided range does not contain the maximum, the returned value
will be at the border.
**/ 
auto find_rhoc_tov_max_mass(eos_barotr eos, 
                       const real_t rhobr0, const real_t rhobr1,
                       const int bits=28, const real_t acc=1e-8, 
                       unsigned int max_steps=30)
-> real_t;


/**\brief Find TOV model with gravitational mass
@param eos The (barotropic) EOS of the NS
@param mg Gravitational mass
@param rhobr0 Lower bound of central density for search
@param rhobr1 Upper bound of central density for search
@param acc Accuracy for TOV solving and root finding
@param max_steps Maximum steps for root finding

@return Central baryonic mass density of TOV model with given mass

This finds a TOV solution of given gravitational mass. It is up to the
user to provide an initial bracket for the central density, such
that the mass as function of density is monotonic (beware the 
different TOV solution branches).
**/ 

auto find_rhoc_tov_of_mass(eos_barotr eos, real_t mg, 
                        const real_t rhobr0, const real_t rhobr1,
                        real_t acc=1e-8, unsigned int max_steps=30)
-> real_t;



}


#endif
