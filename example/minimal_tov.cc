#include "reprimand/unitconv.h"
#include "reprimand/eos_barotr_poly.h"
#include "reprimand/spherical_stars.h"

using namespace EOS_Toolkit;

int main() {
    auto u = units::geom_solar();

    const double n_poly      = 1;
    const double rho_poly    = 6.176e+18 / u.density();
    const double rhomax_poly = 1E40 / u.density();
    auto eos = make_eos_barotr_poly(n_poly, rho_poly,
                                    rhomax_poly);

    const double tov_cnt_rho = 7.9056e+17 / u.density();

    const tov_acc_simple accs{1e-8, 1e-6};

    auto tov = get_tov_star_properties(eos, tov_cnt_rho, accs);

    double gravitational_mass = tov.grav_mass();
    double circumferential_radius = tov.circ_radius();
    double dimless_tidal_deform = tov.deformability().lambda;
}
