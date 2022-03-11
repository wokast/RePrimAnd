
#include<iostream>
#include<string>

#include "unitconv.h"
#include "eos_barotropic.h"
#include "eos_barotr_poly.h"
#include "spherical_stars.h"


using namespace std;
using namespace EOS_Toolkit;


void tovsol_basic(int nrep)
{

  auto u = units::geom_solar();

  const double n_poly      = 1;
  const double rmd_poly    = 6.176e+18 / u.density();
  const real_t rhomax_poly = 1E40 / u.density();
  auto eos = make_eos_barotr_poly(n_poly, rmd_poly, rhomax_poly);

  const double tov_cnt_rmd       = 7.9056e+17 / u.density();
  const tov_acc_simple accs{ 1e-8, 1e-6, 100 };
  
  for (int i=0; i<nrep; ++i) {
    auto tov = make_tov_star(eos, tov_cnt_rmd, accs, false);
  }
}


int main() 
{
  tovsol_basic(10000);  
}


