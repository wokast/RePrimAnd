#ifndef UNITCONV_H
#define UNITCONV_H

#include <iostream>
#include <string>
#include <cmath>

namespace EOS_Toolkit {

///Class to represent physical unit systems
/**
This class describes a system of units, by storing base units for 
time, length and mass, and computes derived units. The unit system 
in which to specify the new units is not fixed, so that this class may 
be used to convert between any two systems of units. The unit object 
can also be used to absolutely specify a unit system. In that case, the 
convention is that the units are expressed in SI units. There are 
predefined unit objects which express different geometric unit 
systems (and SI units) in terms of the SI units.
**/
class units {
  const double ulength{1.0};
  const double utime{1.0};
  const double umass{1.0};
  
  public:
  
  ///Default constructor results in SI units
  constexpr units() = default;
  
  ///Constructor for direct specification of base units
  constexpr units(double ulength_,double utime_,double umass_)
    : ulength(ulength_), utime(utime_), umass(umass_) {}

  ///Express units in terms of other units
  constexpr units operator/(
    const units &base ///<..in terms of those units
  ) const
  {
    return {ulength/base.ulength, utime/base.utime, umass/base.umass}; 
  }
  
  ///Unit of length
  constexpr double length() const {return ulength;}

  ///Unit of time
  constexpr double time() const {return utime;}

  ///Unit of frequency
  constexpr double freq() const {return 1.0/utime;}

  ///Unit of mass
  constexpr double mass() const {return umass;}

  ///Unit of velocity
  constexpr double velocity() const {return ulength/utime;}

  ///Unit of acceleration
  constexpr double accel() const {return velocity()/utime;}

  ///Unit of force
  constexpr double force() const {return accel()*mass();}

  ///Unit of area
  constexpr double area() const {return ulength*ulength;}

  ///Unit of volume
  constexpr double volume() const {return ulength*area();}

  ///Unit of mass density
  constexpr double density() const {return mass()/volume();}

  ///Unit of pressure
  constexpr double pressure() const {return force()/area();}

  ///Unit for moment of inertia
  constexpr double mom_inertia() const {return mass()*area();}

  ///Get SI units
  static constexpr units si() {return {};}

  /**\brief Compute units with G=c=1 and given length unit
  
  @param g_si  Gravitational constant \f$ G \f$ in SI units 
  @return units object
  **/  
  static constexpr units geom_ulength(double ulength, double g_si=G_SI)
  {
    return {ulength, ulength/c_SI, ulength * c_SI * c_SI / g_si};
  }

  /**\brief Compute units with G=c=1 where length unit is meter
  
  @param g_si  Gravitational constant \f$ G \f$ in SI units
  @return units object
  **/
  static constexpr units geom_meter(double g_si=G_SI) 
  {
    return geom_ulength(1.0, g_si);
  }

  /**\brief Compute units with G=c=1 and given density unit
  @param udensity Target density unit in SI units
  @param g_si  Gravitational constant \f$ G \f$ in SI units
  @return units object
  **/  
  static units geom_udensity(double udensity, 
                                       double g_si=G_SI)
  {
    return geom_ulength( c_SI / std::sqrt(g_si * udensity), g_si);
  }

  /**\brief Compute units with G=c=1 and given mass unit

  @param umass Target mass unit in SI units
  @param g_si  Gravitational constant \f$ G \f$ in SI units
  @return units object
  **/  
  static constexpr units geom_umass(double umass, double g_si=G_SI)
  {
    return geom_ulength(umass * g_si /(c_SI*c_SI), g_si);
  }

  /**\brief Compute units with G=c=1 and the mass unit is the solar mass

  @param m_sun_si Solar mass in SI units
  @param g_si  Gravitational constant \f$ G \f$ in SI units
  @return units object
  **/    
  static constexpr units geom_solar(double m_sun_si=M_sun_SI, 
                                    double g_si=G_SI) 
  {
    return geom_umass(m_sun_si, g_si);
  }

  std::string to_str() const;

  ///Speed of light from particle physics data group (2021)
  static constexpr double PDG2021_c_SI = 299792458.0;// m/s

  ///Newtonian graviational constant from particle physics data group (2021)
  static constexpr double PDG2021_G_SI = 6.67430E-11;// m^3 / (kg s^2)
  
  ///Solar mass from particle physics data group (2021)
  static constexpr double PDG2021_M_sun_SI = 1.98841E30;// kg
  
  ///Speed of light
  static constexpr double c_SI = PDG2021_c_SI;
  
  ///Default value used for Newtonian graviational constant
  static constexpr double G_SI = PDG2021_G_SI;
  
  ///Default value used for solar mass
  static constexpr double M_sun_SI = PDG2021_M_sun_SI;
};

std::ostream& operator<<(std::ostream& o, const units& u);

}

#endif

