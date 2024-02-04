#include "reprimand/unitconv.h"
#include "reprimand/eos_barotropic.h"
#include "reprimand/eos_barotr_file.h"
#include "reprimand/star_sequence.h"
#include "reprimand/star_seq_file.h"


using namespace EOS_Toolkit;

int main() {
    // Choose geometric units to work in
    auto u = units::geom_solar();

    //Load EOS from file
    auto eos = load_eos_barotr("example.eos.h5", u);
    
    // Specification for desired accuracy: using defaults
    const auto acc = star_acc_simple();
    
    // Create stable TOV branch
    auto seq = make_tov_branch_stable(eos, acc);
    
    // Save branch to file 
    save_star_branch("example.tovseq.h5", seq);
    
    // Compute tidal deformability at 90% of maximum mass
    
    double mass = 0.9 * seq.grav_mass_maximum();
    
    double lambda09 
      = seq.lambda_tidal_from_grav_mass(mass);
}
