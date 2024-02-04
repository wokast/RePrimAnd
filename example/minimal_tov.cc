#include "reprimand/unitconv.h"
#include "reprimand/eos_barotr_poly.h"
#include "reprimand/spherical_stars.h"

using namespace EOS_Toolkit;

int main() {
    // Choose geometric units to work in
    auto u = units::geom_solar();

    // Create polytropic EOS 
    const double n_poly      = 1;
    const double rho_poly    = 6.176e+18 / u.density();
    const double rhomax_poly = 1E40 / u.density();
    auto eos = make_eos_barotr_poly(n_poly, rho_poly,
                                    rhomax_poly);

    // Choose NS central density
    const double tov_cnt_rho_SI = 7.9056e+17;
    const double tov_cnt_rho = tov_cnt_rho_SI / u.density();

    // Specification for desired accuracy
    const auto accs = star_acc_simple(true,  // want tidal deform Lambda
                                      false, // don't need bulk radius
                                      1e-6,  // accuracy for M,R,I
                                      1e-4   // accuracy for Lambda
                                      );

    // Compute NS properties
    auto tov = get_tov_properties(eos, tov_cnt_rho, accs);

    double gravitational_mass = tov.grav_mass();
    double circumferential_radius = tov.circ_radius();    
    double dimless_tidal_deform = tov.deformability().lambda;
    
    // Convert back units    
    double gravitational_mass_SI = gravitational_mass * u.mass();
    double circumferential_radius_SI = circumferential_radius * u.length();
    
}
