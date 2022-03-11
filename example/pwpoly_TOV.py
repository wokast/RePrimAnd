"""This example demonstrate the use of the Python bindings of the 
RePrimAnd TOV solver and creation of piecewise polytropic EOS on the 
fly (not using EOS files).
"""

import pyreprimand as pyr
import numpy as np
import matplotlib
matplotlib.use("AGG")
from matplotlib import pyplot as plt

def get_seq(eos, rho_min, rho_max,   
            nsamp=400, acc_tov=1.e-6, acc_def=1.e-3, acc_mres=500):
    """
    Compute mass, radius, and tidal deformability for a sequence
    of TOV solution with given EOS and within given range of central
    density.
    """    
    rho_max = min(rho_max, eos.range_rho.max*0.99999)
    
    rhoc = np.exp(np.linspace(np.log(rho_min), np.log(rho_max), nsamp))
    
    accs = pyr.tov_acc_simple(acc_tov, acc_def, acc_mres)
    
    @np.vectorize
    def props(rho_cent):
        tov  = pyr.get_tov_star_properties(eos, rho_cent, accs)
        return (tov.grav_mass, tov.circ_radius, 
                  tov.deformability.lambda_tidal)
        
    return props(rhoc)
    
    
def plot_seqs(seqs):
    """
    Plot mass-radius and compactness-deformability diagrams for 
    a list of TOV sequences.
    """
    u   = pyr.units.geom_solar
    
    fig,(ax1,ax2) = plt.subplots(2,1)
    cm = plt.get_cmap("viridis")
    for oga,mg,rc,lt in seqs:
        clr = cm(oga) 
        c = mg/rc
        ax2.semilogy(c, lt, ls='-', color=clr)
        ax2.set_xlabel(r"$M/R$")
        ax2.set_ylabel(r"$\Lambda$")
        ax1.plot(rc * u.length / 1e3, mg, ls='-', color=clr)
        ax1.set_xlabel(r"$R\,[\mathrm{km}]$")
        ax1.set_ylabel(r"$M \,[M_\odot]$")
    plt.tight_layout()
    plt.savefig("pwpoly_TOV.pdf")
  


def main():
    """
    Set up a family of piecewise polytropic EOS by varying one
    of the segment's adiabatic exponent. Compute TOV sequence
    for each EOS, plot them, save plot in current directory.
    """
    u   = pyr.units.geom_solar
    
    seg_bounds_si = [0.0, 2.44034e+10, 3.78358e+14, 2.6278e+15,
                     9.417030181375906e+16, 5.011872336272714e+17, 
                     1e+18]
    seg_gammas = [1.58425, 1.28733, 0.62223, 
                  1.35692, 3.224, 3.033, 1.325]
    rho_poly_si = 9.535102289214910e+16
    rho_max_si = 2e+18
    rho_min_si = 2e+17
    
    
    seg_bounds = np.array(seg_bounds_si) / u.density
    rho_poly    = rho_poly_si / u.density
    rho_max     = rho_max_si / u.density
    rho_min     = rho_min_si / u.density
    
    
    seqs = []
    for oga in np.linspace(0, 1, 10):
      ga = 2.9 + (oga-0.5)*0.5
      gammas = seg_gammas[:4] + [ga] + seg_gammas[5:]
      eos = pyr.make_eos_barotr_pwpoly(rho_poly, seg_bounds, gammas,
                                     rho_max*2)
                                     
      seq = get_seq(eos, rho_min, rho_max)
      seqs.append([oga]+list(seq))
    
    plot_seqs(seqs)
    

if __name__ == '__main__':
    main()
