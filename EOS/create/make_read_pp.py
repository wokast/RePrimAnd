#!/bin/env python
import numpy as np
import pyreprimand as pyr

C_SI = 299792458.0
CGS_DENS = 1e3
CGS_PRESS = 0.1


def prep_sly():
    

    def matching_ds(rho_bnd, ds_j, ga_j, ga_i):
      n   = 1.0/(ga_i-1.0)
      n_j   = 1.0/(ga_j-1.0)
      et  = n / n_j
      ds  = (ds_j**et) * (rho_bnd**(1.0-et))
      return ds
      
    sly_ga        = np.array([ 
      1.58425,  
      1.28733,  
      0.62223,  
      1.35692
    ])

    sly_rho_b_cgs = np.array([  
      2.44034000e+07,   
      3.78358000e+11,   
      2.62780000e+12
    ])
    sly_rho_b_si = sly_rho_b_cgs * CGS_DENS
    
    sly_k0_cgs = 6.80110e-9
    
    sly_ds_si = [((sly_k0_cgs)**(1/(1-sly_ga[0]))) * CGS_DENS]
    
    for l in range(3):
        ds = matching_ds(sly_rho_b_si[l], sly_ds_si[-1],
                         sly_ga[l],  sly_ga[l+1])
        sly_ds_si.append(ds)
    
    
    return sly_ds_si, sly_ga, sly_rho_b_si



def make_read_eos(name, p1_cgs, ga, rho_max_si=1e20):
    
    p1_si = p1_cgs * CGS_PRESS
    
    sly_ds_si, sly_ga, sly_rho_b_si = prep_sly()
    
    rho_b_si     = np.array([10**14.7, 1e15]) * CGS_DENS
    
    dsx = sly_ds_si[-1]
    gax = sly_ga[-1]
    rho_m_si = dsx * ((dsx * (C_SI**2)/p1_si)
                          * (rho_b_si[0]/dsx)**ga[0] )**(1./(ga[0]-gax))
    
    rho_b_si = np.hstack(([0], sly_rho_b_si, [rho_m_si], rho_b_si))
    ga = np.hstack((sly_ga, ga))

   
    path     = f"{name}.eos.h5"
    
    uc = pyr.units.geom_solar()
    rho_poly = sly_ds_si[0] / uc.density
    rho_b    = rho_b_si / uc.density 
    rho_max = rho_max_si / uc.density
    
    eos = pyr.make_eos_barotr_pwpoly(rho_poly, rho_b, ga, rho_max, uc)
    pyr.save_eos_barotr(path, eos, "A piecewise polytropic EOS from [Read et al]")
    
    rgrho2 = pyr.range(rho_b[1] / 2, eos.range_rho.max / 1.0000001)
    polyn2 =  1./ (ga[0] - 1.)
    ppm2   = 200
    eos2   = pyr.make_eos_barotr_spline(eos, rgrho2, polyn2, ppm2) 
    path2  = f"{name}.spline.eos.h5"
    pyr.save_eos_barotr(path2, eos2, "A spline representation of a piecewise polytropic EOS from [Read et al]")
    
#    

def make_read_cat():
    params = {
      'APR4':(10**34.269, [2.830, 3.445, 3.348]),
      'APR3':(10**34.392, [3.166, 3.573, 3.281]),
      'H4':  (10**34.669, [2.909, 2.246, 2.144]),
      'ALF2':(10**34.616, [4.070, 2.411, 1.890]),
      'MPA1':(10**34.495, [3.446, 3.572, 2.887]),
      'MS1': (10**34.858, [3.224, 3.033, 1.325]),
      'MS1B':(10**34.855, [3.456, 3.011, 1.425]),
      'WFF1':(10**34.031, [2.519, 3.791, 3.660]),
      'WFF2':(10**34.233, [2.888, 3.475, 3.517]),
      'ENG': (10**34.437, [3.514, 3.130, 3.168])
    }

    for n,p in params.items():
        name = f"{n}_Read_PP"
        make_read_eos(name, *p)


make_read_cat()   
