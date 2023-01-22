#ifndef STAR_SEQUENCE_H
#define STAR_SEQUENCE_H
#include <memory>
#include <vector>
#include <cassert>
#include <boost/optional.hpp>
#include "config.h"
#include "intervals.h"
#include "spherical_stars.h"
#include "eos_barotropic.h"
#include "datastore.h"

namespace EOS_Toolkit {
  

namespace detail {
  class star_seq_impl;
  class star_branch_impl;
}  
  
/**\brief Class representing sequences of neutron stars or similar.

This class allows to store precomputed properties for a sequence 
of spherical stars and provide those as functions of the central 
pseudo-enthalpy. For this, regular spaced monotonic spline 
interpolation is used. The unit system can be chosen when creating
sequences, but is assumed to be geometric. It is stored for 
bookkeeping.

Sequence objects can be cheaply copied since they only store a 
reference-counted pointer to the actual (immutable) data. 
Using the objects is thread-save. 
*/   
class star_seq {
  protected:
  
  using impl_t = detail::star_seq_impl;
  using spimpl_t = std::shared_ptr<const impl_t>;
  
  spimpl_t pimpl;
  auto valid() const -> const impl_t&;

  auto implementation() const -> const impl_t& {return valid();}
  
  public:
  
  

  using range_t  = interval<real_t>;
  
  /**\brief Constructor is not intended for direct use
  */
  explicit star_seq(spimpl_t seq) 
  : pimpl(std::move(seq)) 
  {
    assert(pimpl); 
  }
  
  /**\brief Default constructor

  Creates uninitialized object. Any attemt to use it throws an 
  exception. One can use it after assigning an initialized one to it.
  */
  star_seq()                       = default;  
  
  ///Copy constructor.
  star_seq(const star_seq &)       = default;
  
  ///Move constructor
  star_seq(star_seq &&)            = default;
  
  ///Assignment operator.
  star_seq& operator=(const star_seq&)  = default;
  
  ///Move assignment operator
  star_seq& operator=(star_seq&&) = default;  
  
  /**\brief Compute gravitational mass
  
  @param gm1c  Central pseudo enthalpy \f$ g-1 \f$
  @returns    Gravitational mass of the star 
  */
  auto grav_mass_from_center_gm1(real_t gm1c) 
  const -> real_t; 

  /**\brief Compute baryonic mass
  
  @param gm1c  Central pseudo enthalpy \f$ g-1 \f$
  @returns    Baryonic mass of the star 
  */
  auto bary_mass_from_center_gm1(real_t gm1c) const -> real_t; 
  
  /**\brief Compute circumferential radius
  
  @param gm1c  Central pseudo enthalpy \f$ g-1 \f$
  @returns    Proper circumferential radius of the star
  */  
  auto circ_radius_from_center_gm1(real_t gm1c) const -> real_t; 
  
  /**\brief Compute moment of inertia
  
  @param gm1c  Central pseudo enthalpy \f$ g-1 \f$
  @returns    Moment of inertia of the star 
  */
  auto moment_inertia_from_center_gm1(real_t gm1c) const -> real_t; 

  /**\brief Compute dimensionless tidal deformability
  
  @param gm1c  Central pseudo enthalpy \f$ g-1 \f$
  @returns    Dimensionless tidal deformability 
  */
  auto lambda_tidal_from_center_gm1(real_t gm1c) const -> real_t; 
    
  /**\brief Range of central pseudo-enthalpy
  
  @returns The avaialable range of the pseudo-enthalpy \f$ g-1 \f$
  */
  auto range_center_gm1() const -> range_t;
  
  /**\brief Check if given central pseudo enthalpy is within 
  available range
  
  @param gm1c  Central pseudo enthalpy \f$ g-1 \f$
  @returns    If it is covered by the sequence
  */  
  auto contains_gm1(real_t gm1c) const -> bool;

  /**\brief Save sequence to a datastore
  
  Helper function used internally by the sequence file functionality.
  */
  void save(datasink s) const;

  /**\brief Load sequence from a datastore
  
  Helper function used internally by the sequence file functionality.
  */
  auto static from_datasource(datasource s, units u) -> spimpl_t;
  
  /**\brief Obtain unit system of sequence
  
  @returns The (geometric) unit system used for the sequence. 
  */
  auto units_to_SI() const -> const units&;
};
  

/**\brief Class representing a stable branch of a star sequence.

It inherits the interface of star_seq. In addition, it provides
star properties as function of gravitational mass.
*/     
class star_branch : public star_seq {
  private:
  
  using impl_t = detail::star_branch_impl;
  using spimpl_t = std::shared_ptr<const impl_t>;
  
  spimpl_t pimpl;
  const impl_t& valid() const;

  auto implementation() const -> const impl_t& {return valid();}
  
  public:
  using star_seq::range_t;
  
  /**\brief Constructor is not intended for direct use
  */
  explicit star_branch(star_seq::spimpl_t seq, spimpl_t brnch) 
  : star_seq{std::move(seq)}, pimpl{std::move(brnch)} 
  {
    assert(pimpl); 
  }
  
  /**\brief Default constructor

  Creates uninitialized object. Any attemt to use it throws an 
  exception. One can use it after assigning an initialized one to it.
  */
  star_branch()                              = default;  
  
  ///Copy constructor.
  star_branch(const star_branch &)       = default;
  
  ///Move constructor
  star_branch(star_branch &&)            = default;
  
  ///Assignment operator.
  star_branch& operator=(const star_branch&)  = default;
  
  ///Move assignment operator
  star_branch& operator=(star_branch&&) = default;  
  

  /**\brief Compute central pseudo enthalpy from gravitational mass
  
  @param mg   Gravitational mass
  @returns    Central pseudo enthalpy \f$ g-1 \f$ 
  */  
  auto center_gm1_from_grav_mass(real_t mg) const -> real_t; 
  
  /**\brief Compute baryonic mass from gravitational mass
  
  @param mg   Gravitational mass
  @returns    Baryonic mass
  */  
  auto bary_mass_from_grav_mass(real_t mg) const -> real_t; 
  
  /**\brief Compute circumferential radius from gravitational mass
  
  @param mg   Gravitational mass
  @returns    Proper circumferential radius
  */    
  auto circ_radius_from_grav_mass(real_t mg) const -> real_t; 
  
  /**\brief Compute moment of inertia from gravitational mass
  
  @param mg   Gravitational mass
  @returns    Moment of inertia
  */      
  auto moment_inertia_from_grav_mass(real_t mg) const -> real_t; 
  
  /**\brief Compute tidal deformability from gravitational mass
  
  @param mg   Gravitational mass
  @returns    Dimensionless tidal deformability 
  */        
  auto lambda_tidal_from_grav_mass(real_t mg) const -> real_t; 

  /**\brief Range of central pseudo-enthalpy
  
  The upper bound is the central pseudo-enthalpy of the maximum mass
  model unless it was not within the EOS validity range. This
  can be checked with includes_maximum().
  
  @returns The avaialable range of the pseudo-enthalpy \f$ g-1 \f$
  */    
  //Should this be virtual?  
  auto range_center_gm1() const -> range_t;
  
  /**\brief Check if given central pseudo enthalpy is within 
  available range
  
  @param gm1c  Central pseudo enthalpy \f$ g-1 \f$
  @returns    If it is covered by the sequence
  */    
  auto contains_gm1(real_t gm1c) const -> bool;

  /**\brief Compute gravitational mass from central pseudo enthalpy
  
  @param gm1c  Central pseudo enthalpy \f$ g-1 \f$
  @returns    Gravitational mass of the star 
  */
  auto grav_mass_from_center_gm1(real_t gm1c) const -> real_t; 
  
  /**\brief Compute baryonic mass from central pseudo enthalpy
  
  @param gm1c  Central pseudo enthalpy \f$ g-1 \f$
  @returns    Baryonic mass of the star 
  */
  auto bary_mass_from_center_gm1(real_t gm1c) const -> real_t; 
  
  /**\brief Compute circumferential radius from central pseudo enthalpy
  
  @param gm1c  Central pseudo enthalpy \f$ g-1 \f$
  @returns    Proper circumferential radius of the star
  */
  auto circ_radius_from_center_gm1(real_t gm1c) const -> real_t; 
  
  /**\brief Compute moment of inertia from central pseudo enthalpy
  
  @param gm1c  Central pseudo enthalpy \f$ g-1 \f$
  @returns    Moment of inertia of the star 
  */
  auto moment_inertia_from_center_gm1(real_t gm1c) const -> real_t; 
  
  /**\brief Compute dimensionless tidal deformability 
  from central pseudo enthalpy
  
  @param gm1c  Central pseudo enthalpy \f$ g-1 \f$
  @returns    Dimensionless tidal deformability 
  */
  auto lambda_tidal_from_center_gm1(real_t gm1c) const -> real_t; 
  
  /**\brief Available range of gravitational mass
  
  The upper bound is the maximum mass unless the maximum mass model
  central density was not within the EOS validity range. This
  can be checked with includes_maximum().
  
  @returns   Gravitational mass range
  */    
  auto range_grav_mass() const -> range_t;
  
  /**\brief Check if given gravitational mass is within 
  available range
  
  @param mg   Gravitational mass
  @returns    If it is covered by the sequence
  */  
  auto contains_grav_mass(real_t mg) const -> bool;

  
  auto includes_maximum() const -> bool;
  auto grav_mass_maximum() const -> real_t;
  auto bary_mass_maximum() const -> real_t;
  auto center_gm1_maximum() const -> real_t;

  /**\brief Save branch to a datastore
  
  Helper function used internally by the sequence file functionality.
  **/
  void save(datasink s) const;
  
  /**\brief Load branch from a datastore
  
  Helper function used internally by the sequence file functionality.
  **/
  auto static from_datasource(datasource s, units u) -> spimpl_t;

};

/**\brief Create sequence of NS properties from existing data

@param mg Sample points of gravitational mass
@param mb Sample points of baryonic mass
@param rc Sample points of proper circumferential radius
@param mi Sample points of moment of inertia
@param lt Sample points of dimensionless tidal deformability
@param rg_gm2 Range of central pseudo enthalpy
@param u Unit system of the sequence, assumed geometric.

@return A star_seq object describing a NS sequence. 

for this function, the sample points must be uniformly spaced in 
the central pseudo enthalpy.
**/ 
auto make_star_seq(std::vector<real_t> mg, 
       std::vector<real_t> mb, std::vector<real_t> rc, 
       std::vector<real_t> mi, std::vector<real_t> lt, 
       star_seq::range_t rg_gm1, units u)
-> star_seq;


/**\brief Compute sequence of TOV solutions

@param eos The (barotropic) EOS of the NSs. 
@param acc Tolerances for adaptive ODE solver.
@param rg_gm1 Range of central pseudo enthalpy \f$ g-1 \f$
@param num_samp Sample resolution of the sequence.

@return A star_seq object describing the TOV sequence. 
**/ 
auto make_tov_seq(eos_barotr eos, tov_acc_simple acc, 
                       interval<real_t> rg_gm1, 
                       unsigned int num_samp=500)
-> star_seq;


/**\brief Compute stable branch of TOV solutions

@param eos The (barotropic) EOS of the NSs. 
@param acc Tolerances for adaptive ODE solver.
@param mgrav_min Minimum mass which sequence should cover.
@param num_samp Sample resolution of the sequence.
@param gm1_initial Central enthalpy to indicate desired branch.
@param max_margin Defines when maximum is considerd physical.

@return A star_branch onject describing the TOV stable branch. 

This function employs heuristic algorithm to find the stable branch
of TOV solutions. Since there may be more than one such branch, one 
has to provide a central pseudo-enthalpy to indicate the correct one.
By increasing/decreasing this value successively by some factor,
a search interval is expanded until it brackets the maximum mass
or until the upper bound exceeds the EOS validity range. The initial 
guess is not required to be within a stable branch, and the default
value should work for any remotely realistic NS EOS.
In a similar fashion, a central pseudo-enthalpy for which the NS mass 
falls below the parameter mgrav_min is determined. Next, the maximum
is determined using a maximum search. However, the maximum might be 
located at the EOS validity bound. The maxmimum is considered physical 
based on a simple heuristics: (g_max-1)*(1+max_margin) < g_eos, where 
g_max and g_eos are the the central pseudo-enthalpy of the maximum mass 
model and the EOS upper validity bound. This criterion can later be 
queried using the includes_maximum() method of the branch object.
Finally, the branch is then sampled with resolution given by num_samp.
**/ 
auto make_tov_branch_stable(eos_barotr eos, tov_acc_simple acc,
                          real_t mgrav_min=0.5,
                          unsigned int num_samp=500,
                          real_t gm1_initial=1.2, 
                          real_t max_margin=5e-2)
-> star_branch;

}

#endif

