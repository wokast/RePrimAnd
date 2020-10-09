#include <boost/test/unit_test.hpp>

#include "test_utils.h"
#include <algorithm>
#include <boost/format.hpp>

using boost::format;
using boost::str;


failcount::~failcount()
{
  BOOST_CHECK_MESSAGE(count==0, 
                 count << " fails while testing: " << descr );
}

void failcount::fail(std::string msg) 
{
  ++count;
  BOOST_CHECK_MESSAGE(false, "Violated expectation: " << msg);
}

bool failcount::istrue(bool c, std::string msg) 
{
  if (!c) fail(msg);
  return c;
}

bool failcount::isnan(double v, std::string msg) 
{
  return istrue(std::isnan(v), msg + " is NAN");
}

bool failcount::isfinite(double v, std::string msg) 
{
  return istrue(std::isfinite(v), msg + " is finite");
}

bool failcount::isless(double smaller, double larger, std::string msg)
{
  return istrue(smaller < larger, msg); 
}

bool failcount::isleq(double smaller, double larger, std::string msg)
{
  return istrue(smaller <= larger, msg); 
}

bool failcount::isclose(double a, double b, 
                        double reltol, double abstol, std::string msg)
{
  bool c= istrue(isfinite(reltol, "rel. tol.") 
                    && isfinite(abstol, "abs. tol."),
                  msg + " comparison tolerances finite");
  c= c && istrue(isfinite(a, "first") && isfinite(b, "second"), 
                  msg + " compared values finite");
  if (c) {
    c = istrue(fabs(a-b) <= std::max(fabs(a),fabs(b)) * reltol + abstol, 
                msg + " within tolerance"); 
    if (!c) {
      double e = (a-b) * 2.0/(fabs(a)+fabs(b));
      postmortem(str(format("a=%.15e, b=%.15e, err = %.3g") 
                      % a % b % e)); 
    }
  }

  return c;
}                        

void failcount::postmortem(std::string msg) const
{
  BOOST_CHECK_MESSAGE(false, msg);
}

