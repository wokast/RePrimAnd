#include <iostream>
#include "reprimand/con2prim_imhd.h"
#include "reprimand/eos_idealgas.h"

using namespace EOS_Toolkit;


int main()
{
  //Get some EOS
  real_t max_eps = 11.;
  real_t max_rho = 1e6;
  real_t adiab_ind = 1.0;
  auto eos = make_eos_idealgas(adiab_ind, max_eps, max_rho);
  
  //Set up atmosphere
  real_t atmo_rho = 1e-20;
  real_t atmo_eps = 0.1;
  real_t atmo_ye = 0.5;
  real_t atmo_cut = atmo_rho * 1.01;
  real_t atmo_p = eos.at_rho_eps_ye(atmo_rho, atmo_eps, atmo_ye).press();

  atmosphere atmo{atmo_rho, atmo_eps, atmo_ye, atmo_p, atmo_cut};

  //Primitive recovery parameters 
  real_t rho_strict = 1e-11;
  bool  ye_lenient = false;
  int max_iter = 30;
  real_t c2p_acc = 1e-8;
  real_t max_b = 10.;
  real_t max_z = 1e3;
  
  //Get a recovery function
  con2prim_mhd cv2pv(eos, rho_strict, ye_lenient, max_z, max_b, 
                     atmo, c2p_acc, max_iter);


  //Some example values
  real_t dens = 1e-6;
  real_t tau  = 1e-8;
  real_t scon_x = dens * 5.;
  real_t bcon_z = 1e-4 * sqrt(dens);
  real_t trye = 0.5*dens;


  //Imagine this part loops over your evolution grid
  {     
    //collect
    cons_vars_mhd evolved{dens, tau, trye, 
                          {scon_x,0.,0.}, {0.,0.,bcon_z}};    
    sm_metric3 g;
    g.minkowski();

    prim_vars_mhd primitives;
    con2prim_mhd::report rep;
    
    //recover
    cv2pv(primitives, evolved, g, rep);
    
    //check
    if (rep.failed())
    {
      std::cerr << rep.debug_message(); 
      //abort simulation
    }
    else {
      //write back primitive vars to grid here
    
      if (rep.adjust_cons) {
        //write back corrected evolved vars to grid here
      }
    }
  }
  
  return 0;
}
