"""This script is run to test that the Python bindings of the 
RePrimAnd TOV solver C++ API can be imported and executed."""

import pyreprimand as pyr
import numpy as np
import sys
import pathlib

def allgood(path, eosname="H4_Read_PP"):
    fn = path/f"{eosname}.eos.h5"
    try:
        u   = pyr.units.geom_solar

        cnt_rho = 7.9056e+17 / u.density
        
        eos = pyr.load_eos_barotr(str(fn), units=u)
        
        accs = pyr.tov_acc_simple(1e-8, 1e-6)
        tov1 = pyr.get_tov_star_properties(eos, cnt_rho, accs)
        
        mg = tov1.grav_mass
        mb = tov1.bary_mass
        rc = tov1.circ_radius
        pv = tov1.proper_volume
        mi = tov1.moment_inertia
        lt = tov1.deformability.lambda_tidal
        k2 = tov1.deformability.k2
        
        tov2 = pyr.make_tov_star(eos, cnt_rho, accs)
        
        mg = tov2.grav_mass
        mb = tov2.bary_mass
        rc = tov2.circ_radius
        pv = tov2.proper_volume
        mi = tov2.moment_inertia
        lt = tov2.deformability.lambda_tidal
        k2 = tov2.deformability.k2
        
        
        r = np.linspace(0,30e3/u.length, 100)
        
        rho  = tov2.rho_from_rc(r)
        eps  = tov2.eps_from_rc(r)
        rhoE = rho * (1+eps)
        p    = tov2.press_from_rc(r)
        cs   = tov2.csnd_from_rc(r)
        
    
      
    except Exception as e:
        print(e)
        return False
        
    return True
#


sys.exit(0 if allgood(pathlib.Path(sys.argv[1])) else 1)




