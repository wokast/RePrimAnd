"""This script is run to test that the Python bindings of the 
RePrimAnd Barotropic EOS C++ API can be imported and executed."""

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
        
        rho = np.linspace(0,eos.range_rho.max, 100)
        
        eps = eos.eps_at_rho(rho)
        p   = eos.press_at_rho(rho)
        cs  = eos.csnd_at_rho(rho)
        gm1 = eos.gm1_at_rho(rho)
        
        rho = eos.rho_at_gm1(gm1)
        eps = eos.eps_at_gm1(gm1)
        p   = eos.press_at_gm1(gm1)
        cs  = eos.csnd_at_gm1(gm1)
    
      
    except Exception as e:
        print(e)
        return False
        
    return True
#


sys.exit(0 if allgood(pathlib.Path(sys.argv[1])) else 1)




