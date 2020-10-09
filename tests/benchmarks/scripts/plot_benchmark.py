#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from stdplt import *

GAUSS_GU = 1.1973228161170991e-20

def plot_perf_zeps(path, eosname, dat_bz, dat_bl):
  figname = "perf_eos%s_z_eps" % eosname
  
  with Figure(os.path.join(path, figname)) as fig:
    grid = twopanel(fig)
    
    maxit     = max(int(dat_bz[-1][2].max()), int(dat_bl[-1][2].max()))
    bticks    = np.arange(1,maxit, max(1,int(maxit/10)))
    clrmap    = get_color_map("viridis", over='r')

    for (x0,x1,(z,eps,it,p,B)),ax in zip([dat_bz, dat_bl], grid):
      im = plot_color(it, x0, x1, vmin=0.999, vmax=maxit,
                 cmap=clrmap, bar=False, axes=ax)
      ax.set_aspect('auto')
      ax.set_ylabel(r'$\log_{10}(\epsilon_\mathrm{th})$')
    grid[1].set_xlabel(r'$\log_{10}(z[c])$')
    color_bar(im, cax=ax.cax, ticks=bticks, label='Calls to EOS', 
              extend='max')
  #
#

def plot_perf_zb(path, eosname, dat_c, dat_h):
  figname = "perf_eos%s_z_b" % eosname
  with Figure(os.path.join(path, figname)) as fig:
    grid = twopanel(fig)
    
    maxit    = max(int(dat_c[-1][2].max()), int(dat_h[-1][2].max()))
    bticks   = np.arange(1,maxit, max(1,int(maxit/10)))
    clrmap    = get_color_map("viridis", over='r')
    for (x0,x1,(z,b,it,p,B)),ax in zip([dat_c, dat_h], grid):
      im = plot_color(it, x0, x1, vmin=0.999, vmax=maxit, 
                      cmap=clrmap, bar=False, axes=ax)
      
      z1d = z[:,0]
      b1d = b[0,:]
      lgbeta = np.log10(B**2 / p)
      cs = ax.contour(np.log10(z1d), np.log10(b1d), lgbeta.T, 
                   levels=np.arange(-2,11,4), colors='y')

      ax.clabel(cs, inline=True)
      ax.contour(np.log10(z1d), np.log10(b1d), B.T, 
                   levels=[1e16*GAUSS_GU], colors=['r'])

      ax.set_aspect('auto')
      ax.set_ylabel(r'$\log_{10}(B / \sqrt{D})$')
    grid[1].set_xlabel(r'$\log_{10}(z[c])$')
    color_bar(im, cax=ax.cax, ticks=bticks, label='Calls to EOS', 
              extend='max')
  #
#

def load_data(path, eos):
    def read(t):
      c0,c1,d = load_grid(os.path.join(path, t % eos), 5)
      return c0,c1,d
    #
    files = ["perf_eos%s_z_b_cold.dat", "perf_eos%s_z_b_hot.dat",
             "perf_eos%s_z_eps_Bzero.dat", 
             "perf_eos%s_z_eps_Blarge.dat"]
    return map(read, files)
#


def main():
  path = sys.argv[1]
  for eos in ['ig', 'hyb']:
    dat_c, dat_h, dat_bz, dat_bl = load_data(path, eos)
    plot_perf_zeps(path, eos, dat_bz, dat_bl)
    plot_perf_zb(path, eos, dat_c, dat_h)
  #
#

main()
