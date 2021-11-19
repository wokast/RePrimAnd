#include "unitconv.h"
#include <boost/format.hpp>

using namespace std;
using namespace EOS_Toolkit;


constexpr double units::PDG2021_c_SI;
constexpr double units::PDG2021_G_SI;
constexpr double units::PDG2021_M_sun_SI;
constexpr double units::c_SI;
constexpr double units::G_SI;
constexpr double units::M_sun_SI;

string units::to_str() const
{
  boost::format fmt("ulength=%.15e m, utime=%.15e s, umass=%.15e kg");
  fmt % ulength % utime % umass;
  return fmt.str();
}

ostream& EOS_Toolkit::operator<<(ostream& o, const units& u) 
{
  o << u.to_str();
  return o;
}


