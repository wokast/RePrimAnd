#include "reprimand/unitconv.h"
#include "reprimand/eos_barotropic.h"
#include "reprimand/eos_barotr_file.h"
#include "reprimand/star_sequence.h"


using namespace EOS_Toolkit;

int main() {
    auto u = units::geom_solar();

    auto eos = load_eos_barotr("example.eos.h5", u);
    
    const tov_acc_simple acc{1e-8, 1e-6};
    auto seq = make_tov_branch_stable(eos, acc);
    
    
    double max_mass = seq.grav_mass_maximum();
    
    double lambda09 
      = seq.lambda_tidal_from_grav_mass(0.9*max_mass);
}
