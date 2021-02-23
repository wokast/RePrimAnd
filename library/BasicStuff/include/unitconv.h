#ifndef UNITCONV_H
#define UNITCONV_H

#include <iostream>
#include <string>

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
  double ulength{1.0};
  double utime{1.0};
  double umass{1.0};
  
  public:
  
  units() = default;
  ///Constructor for direct specification of base units
  units(double ulength_,double utime_,double umass_)
    : ulength(ulength_), utime(utime_), umass(umass_) {}

  ///Express units in terms of other units
  units operator/(
    const units &base ///<..in terms of those units
  ) const;

  ///Unit of length
  double length() const {return ulength;}

  ///Unit of time
  double time() const {return utime;}

  ///Unit of frequency
  double freq() const {return 1.0/utime;}

  ///Unit of mass
  double mass() const {return umass;}

  ///Unit of velocity
  double velocity() const {return ulength/utime;}

  ///Unit of acceleration
  double accel() const {return velocity()/utime;}

  ///Unit of force
  double force() const {return accel()*mass();}

  ///Unit of area
  double area() const {return ulength*ulength;}

  ///Unit of volume
  double volume() const {return ulength*area();}

  ///Unit of mass density
  double density() const {return mass()/volume();}

  ///Unit of pressure
  double pressure() const {return force()/area();}

  ///Unit for moment of inertia
  double mom_inertia() const {return mass()*area();}

  ///Get SI units
  static units si() {return {};}

  ///Compute units with G=c=1 and given length unit
  static units geom_ulength(double ulength);

  ///Compute units with G=c=1 where length unit is meter
  static units geom_meter() {return geom_ulength(1.0);}

  ///Compute units with G=c=1 and given density unit
  static units geom_udensity(double udensity);

  ///Compute units with G=c=1 and given mass unit
  static units geom_umass(double umass);

  ///Compute units with G=c=1 and the mass unit is the solar mass
  static units geom_solar() {return geom_umass(M_sun_SI);}

  std::string to_str() const;

  static const double c_SI, G_SI, M_sun_SI;
};

std::ostream& operator<<(std::ostream& o, const units& u);

}

#endif

