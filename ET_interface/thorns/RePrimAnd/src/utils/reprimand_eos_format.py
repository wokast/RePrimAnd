import numpy as np
import h5py

FORMAL_M_BARYON_SI = 1.66E-27 #kg

EOS_FILE_FORMAT_VERSION = 0.1


def eosfile_common_(f5, name, comment, eostype):
    """Write fields common to all EOS types in HDF5 file"""
    f5.attrs['format_version'] = EOS_FILE_FORMAT_VERSION
    f5.attrs['eos_type']       = eostype
    f5.attrs['eos_name']       = name
    f5.attrs['eos_comment']    = comment
#

def write_barotr_table_(f, name, comment, isentropic, poly_n,
         rho, eps, press, csnd, gm1,
         ye=None, temp=None):
    """Write tabulated barotropic EOS data into hdf5 file"""
    
    eos_type = 'barotr_table'
    eosfile_common_(f, name, comment, eos_type)
    g   = f.create_group('eos_%s' % eos_type)
        
    def dset(name, data, descr):
      s = g.create_dataset(name, data=np.array(data).astype(np.float))
      s.attrs['description'] = descr
    #
    
    dset('rmd', rho, "Rest mass density [kg m^-3]")
    dset('sed', eps, "Specific energy [dimensionless]")
    dset('press', press, "Pressure [Pa]")
    dset('gm1', gm1, "Pseudo-enthalpy g-1 [dimensionless]")
    dset('csnd', csnd, "Soundspeed [m/s]")
    
    if ye:
      dset('efr', ye, "Electron fraction [dimensionless]")
    #
    if temp:
      dset('temp', temp, "Temperature [MeV]")
    #
    g.attrs['isentropic'] = int(isentropic)
    g.attrs['mbar_si']    = float(FORMAL_M_BARYON_SI)
    g.attrs['poly_n']     = float(poly_n)
#


def comment_pwpoly_(name, rho_poly, rho_max,
                    rho_bnd, gammas, comment=None):
    """Autogenerate description text for piecewise polytropic EOS"""
    sbd = [("%.15e" % b) for b in rho_bnd]
    sbd = "[%s]" % (', '.join(sbd))
    
    sgs = [("%.15e" % g) for g in gammas]
    sgs = "[%s]" % (', '.join(sgs))
    
    s = """Piecewise polytropic EOS %s
  Allowed range (0, %.15e) kg/m^3 
  Polytropic density scale %.15e kg/m^3
  Segment boundaries %s kg/m^3 
  Segment gammas %s 
  
    """ % (name, rho_max, rho_poly, sbd, sgs)
    
    if comment: 
        s += comment
       
    return s
#
  
def write_barotr_pwpoly_(f, name, comment, rho_poly, rho_max,
                         rho_bnd, gammas):
    """Write piecewise polytropic EOS data into hdf5 file"""
    eos_type = 'barotr_pwpoly'
    
    eosfile_common_(f, name, comment, eos_type)
    g   = f.create_group('eos_%s' % eos_type)

    def dset(name, data, descr):
      s = g.create_dataset(name, data=np.array(data).astype(np.float))
      s.attrs['description'] = descr
    #                      
    
    g.attrs['rho_poly']     = float(rho_poly)
    g.attrs['rho_max']      = float(rho_max)
    
    dset('rho_bound', rho_bnd, 
         "Baryon mass density segment boundaries [kg/m^3]")                         
    dset('gamma', gammas, 
         "Adiabatic exponent gamma of segments [dimensionless]")                         
#

def comment_poly_(name, rho_poly, rho_max, n_adiab, comment=None):
    """Autogenerate description text for polytropic EOS"""
    gamma = 1.0 + 1.0/n_adiab
    
    s = """Polytropic EOS %s
  Allowed range (0, %.15e) kg/m^3 
  Polytropic density scale %.15e kg/m^3
  Gamma: %.15f 
  
    """ % (name, rho_max, rho_poly, gamma)
    
    if comment: 
        s += comment
       
    return s
#
  
def write_barotr_poly_(f, name, comment, rho_poly, rho_max, n_adiab):
    """Write polytropic EOS data into hdf5 file"""
    eos_type = 'barotr_poly'
    
    eosfile_common_(f, name, comment, eos_type)
    g   = f.create_group('eos_%s' % eos_type)
    
    g.attrs['rho_poly']     = float(rho_poly)
    g.attrs['rho_max']      = float(rho_max)
    g.attrs['poly_n']       = float(n_adiab)
#


def common_thermal_hybrid_(f, name, cold_name, cold_comment,
                           gamma_th, eps_max):
    """Create hdf5 file with hybrid thermal EOS"""

    comment = """Hybrid EOS 
Gamma_thermal = %.5f
max_eps       = %.5f
Based on cold EOS %s:

%s
""" % (gamma_th, eps_max, cold_name, cold_comment)
    
    eos_type = 'thermal_hybrid'
    eosfile_common_(f, name, comment, eos_type)
    g  = f.create_group('eos_%s' % eos_type)
    g.attrs['gamma_th'] = float(gamma_th)
    g.attrs['eps_max']  = float(eps_max)
    gc  = g.create_group('eos_cold')
    
    return gc
#

def comment_igas_(name, rho_max, n_adiab, comment):

    gamma = 1. + 1. / n_adiab
    s = """Classical ideal gas EOS %s
  Allowed range (0, %.15e) kg/m^3 
  Gamma: %.15f 
  
    """ % (name, rho_max, gamma)
    
    if comment: 
        s += comment
    
    return s
#

def save_barotr_table(path, name, comment, isentropic, poly_n,
         rho, eps, press, csnd, gm1,
         ye=None, temp=None):
    r"""Create hdf5 file with tabulated barotropic EOS

    :param path:  Path of the file to create. Ending should be .eos.h5
    :param name:  Name of the EOS 
    :param comment: Description of the EOS
    :param isentropic: Whether the EOS is isentropic.    
    :param poly_n: Adiabatic index used for extrapolating below
                   lowest tabulated density
    :param rho: Array with baryonic mass density 
                :math:`\rho\, [\mathrm{kg}/\mathrm{m}^3]` 
                at the sample points
    :param eps: Array with specific internal energy 
                :math:`\epsilon` [dimensionless]
                at sample points
    :param press: Array with pressure :math:`P\,[\mathrm{Pa}]` 
                  at sample points
    :param csnd: Array with adiabatic sound speed 
                 :math:`c_s\,[\mathrm{m}/\mathrm{s}]`
                 at sample points
    :param gm1: Array with pseudo-enthalpy 
                :math:`g-1` [dimensionless]
                at sample points
    :param ye: Optionally, array with electron fraction 
               :math:`Y_e` [dimensionless]
               at sample points.
    :param ye: Optionally, array with temperature 
               :math:`T \,[\mathrm{MeV}]` 
               at sample points.
    """
    with h5py.File(path, mode = "w") as f:    
        write_barotr_table_(f, name, comment, isentropic, poly_n,
         rho, eps, press, csnd, gm1,
         ye=ye, temp=temp)
        
#



def save_barotr_poly(path, name, rho_poly, rho_max, n_adiab, 
                       comment=None):                         
    r"""Create hdf5 file with polytropic barotropic EOS

    :param path:  Path of the file to create. Ending should be .eos.h5
    :param name:  Name of the EOS 
    :param rho_poly: Polytropic density scale 
                     :math:`\rho_p\, [\mathrm{kg}/\mathrm{m}^3]` 
    :param rho_max: Maximum allowed baryonic mass density
            :math:`\rho_\mathrm{max} \, [\mathrm{kg}/\mathrm{m}^3]`
    :param n_adiab: Adiabatic index :math:`n`
    :param comment: Optional additional description of the EOS,  
                    appended to autogenerated description
    """
    cmt = comment_poly_(name, rho_poly, rho_max, n_adiab, comment)
    
    with h5py.File(path, mode = "w") as f:    
        write_barotr_poly_(f, name, cmt, rho_poly, rho_max, n_adiab)
    #
#

def save_barotr_pwpoly(path, name, rho_poly, rho_max,
                       rho_bnd, gammas, comment=None):                         
    r"""Create hdf5 file with piecewise poly barotropic EOS

    :param path:  Path of the file to create. Ending should be .eos.h5
    :param name:  Name of the EOS 
    :param rho_poly: Polytropic density scale 
                     :math:`\rho_p\, [\mathrm{kg}/\mathrm{m}^3]` 
                     of first segment
    :param rho_max: Maximum allowed baryonic mass density
            :math:`\rho_\mathrm{max} \, [\mathrm{kg}/\mathrm{m}^3]`
    :param rho_bnd: Array with lower boundary 
                    :math:`\rho_{b,i}\, [\mathrm{kg}/\mathrm{m}^3]` 
                    of baryonic mass densities for each segment. 
                    Must start at zero.
    :param gammas: Array with adiabatic exponents :math:`\Gamma_i`
                   for each segment
    :param comment: Optional additional description of the EOS,  
                    appended to autogenerated description
    """
    cmt = comment_pwpoly_(name, rho_poly, rho_max, rho_bnd, gammas, 
                          comment)
    
    with h5py.File(path, mode = "w") as f:    
        write_barotr_pwpoly_(f, name, cmt, rho_poly, rho_max,
                             rho_bnd, gammas)
    #
#



def save_thermal_hybrid_table(path, name, cold_name, cold_comment, 
         gamma_th, eps_max, poly_n, rho, eps, press, 
         csnd, gm1, ye=None):
    r"""Create a file for a hybrid EOS based on a tabulated cold EOS
    
    :param path:  Path of the file to create. Ending should be .eos.h5
    :param name:  Name of the EOS 
    :param cold_name:  Name of the underlying cold barotropic EOS
    :param cold_comment: Description of underlying cold EOS    
    :param gamma_th: :math:`\Gamma_{th}` defining the thermal part
    :param eps_max: Maximum allowed specific internal energy 
                    :math:`\epsilon` 
    :param poly_n: Adiabatic index used for extrapolating below
                   lowest tabulated density
    :param rho: Array with baryonic mass density 
                :math:`\rho\, [\mathrm{kg}/\mathrm{m}^3]` 
                at the sample points
    :param eps: Array with specific internal energy 
                :math:`\epsilon` [dimensionless]
                at sample points
    :param press: Array with pressure :math:`P\,[\mathrm{Pa}]` 
                  at sample points
    :param csnd: Array with adiabatic sound speed 
                 :math:`c_s\,[\mathrm{m}/\mathrm{s}]`
                 at sample points
    :param gm1: Array with pseudo-enthalpy 
                :math:`g-1` [dimensionless]
                at sample points
    :param ye: Optionally, array with electron fraction 
               :math:`Y_e` [dimensionless]
               at sample points.
    """
    with h5py.File(path, mode = "w") as f:    
        g = common_thermal_hybrid_(f, name, cold_name, cold_comment, 
                                   gamma_th, eps_max)

        write_barotr_table_(g, cold_name, cold_comment, 
                True, poly_n, rho, eps, press, csnd, gm1, ye=ye)        
   #
#

def save_thermal_hybrid_pwpoly(path, name, cold_name, 
                               gamma_th, eps_max,
                               rho_poly, rho_max, rho_bnd, 
                               gammas, cold_comment=None):
    r"""Create a file for a hybrid EOS based on a piecewise polytropic 
    cold EOS
    
    :param path:  Path of the file to create. Ending should be .eos.h5
    :param name:  Name of the EOS 
    :param cold_name:  Name of the underlying cold barotropic EOS
    :param gamma_th: :math:`\Gamma_{th}` defining the thermal part
    :param eps_max: Maximum allowed specific internal energy 
                    :math:`\epsilon` 
    :param rho_poly: Polytropic density scale 
                     :math:`\rho_p\, [\mathrm{kg}/\mathrm{m}^3]` 
                     of first segment
    :param rho_max: Maximum allowed baryonic mass density
            :math:`\rho_\mathrm{max} \, [\mathrm{kg}/\mathrm{m}^3]`
    :param rho_bnd: Array with lower boundary 
                    :math:`\rho_{b,i}\, [\mathrm{kg}/\mathrm{m}^3]` 
                    of baryonic mass densities for each segment. 
                    Must start at zero.
    :param gammas: Array with adiabatic exponents :math:`\Gamma_i`
                   for each segment
    :param cold_comment: Optionally, additional description of the 
                         cold EOS.
    """
    cmt = comment_pwpoly_(name, rho_poly, rho_max, rho_bnd, gammas, 
                          cold_comment)

    with h5py.File(path, mode = "w") as f:    
        g = common_thermal_hybrid_(f, name, cold_name, cmt, 
                                   gamma_th, eps_max)       
        
        write_barotr_pwpoly_(g, cold_name, cmt, rho_poly, rho_max,
                             rho_bnd, gammas) 
    #

#


def save_thermal_hybrid_poly(path, name, cold_name, 
                    gamma_th, eps_max,
                    rho_poly, rho_max, n_adiab,  cold_comment=None):
    r"""Create a file for a hybrid EOS based on a polytropic cold EOS
    
    :param path:  Path of the file to create. Ending should be .eos.h5
    :param name:  Name of the EOS 
    :param cold_name:  Name of the underlying cold barotropic EOS
    :param gamma_th: :math:`\Gamma_{th}` defining the thermal part
    :param eps_max: Maximum allowed specific internal energy 
                    :math:`\epsilon` 
    :param rho_poly: Polytropic density scale 
                     :math:`\rho_p\, [\mathrm{kg}/\mathrm{m}^3]` 
    :param rho_max: Maximum allowed baryonic mass density
            :math:`\rho_\mathrm{max} \, [\mathrm{kg}/\mathrm{m}^3]`
    :param n_adiab: Adiabatic index :math:`n`
    :param cold_comment: Optionally, additional description of the 
                         cold EOS.
    """
    cmt = comment_poly_(name, rho_poly, rho_max, n_adiab, cold_comment)

    with h5py.File(path, mode = "w") as f:    
        g = common_thermal_hybrid_(f, name, cold_name, cmt, 
                                   gamma_th, eps_max)       
       
        write_barotr_poly_(g, cold_name, cmt, rho_poly, rho_max, 
                           n_adiab) 
    #
#





def save_thermal_idealgas(path, name, n_adiab, eps_max,
                          rho_max, comment=None):
    r"""Create a file for a classical ideal gas EOS 
    
    :param path:  Path of the file to create. Ending should be .eos.h5
    :param name:  Name of the EOS 
    :param n_adiab: Adiabatic index :math:`n`
    :param eps_max: Maximum allowed specific internal energy 
                    :math:`\epsilon` 
    :param rho_max: Maximum allowed baryonic mass density
            :math:`\rho_\mathrm{max} \, [\mathrm{kg}/\mathrm{m}^3]`
    :param comment: Optional additional description of the EOS,  
                    appended to autogenerated description
    """
    cmt = comment_igas_(name, rho_max, n_adiab, comment)
    
    with h5py.File(path, mode = "w") as f:    
        eos_type = 'thermal_idealgas'
        eosfile_common_(f, name, cmt, eos_type)
        g   = f.create_group('eos_%s' % eos_type)
    
        g.attrs['eps_max']      = float(eps_max)
        g.attrs['rho_max']      = float(rho_max)
        g.attrs['adiab_index']  = float(n_adiab)    
    #
#






