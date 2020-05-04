#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
from stdplt import *




def plot_acc_zeps(path, eosname, dat_bz, dat_bl):
  figname = "acc_eos%s_z_eps" % eosname
  dats = [dat_bz, dat_bl]
  with Figure(os.path.join(path, figname)) as fig:
    grid = twopanel(fig)
    clrmap    = get_color_map("viridis", over='r', under='k')
    
    maxerr = max((np.max(d[8]) for x0,x1,d in dats)) * 1.001
    for (x0,x1,dat),ax in zip(dats, grid):
      err = np.maximum(dat[8], 1e-20) 
      
      im = plot_color(np.log10(err), x0, x1, 
                      vmin=-16, vmax=log10(maxerr),
                      cmap=clrmap, bar=False, axes=ax)
      
      rnd = dat[10]
      z1d = dat[0][:,0]
      th1d = dat[1][0,:]
      lvl = 1.0
      if np.max(rnd) > lvl:
        ax.contour(np.log10(z1d), np.log10(th1d), rnd.T, 
                   [lvl], colors=['k'])

      ax.set_aspect('auto')
      ax.set_ylabel(r'$\log_{10}(\epsilon_\mathrm{th})$')
    grid[1].set_xlabel(r'$\log_{10}(z[c])$')
    color_bar(im, cax=ax.cax, label=r'$\log_{10}$(relative error)', extend='both')
  #
#

def plot_acc_zb(path, eosname, dat_c, dat_h):
  figname = "acc_eos%s_z_b" % eosname
  dats = [dat_c, dat_h]
  with Figure(os.path.join(path, figname)) as fig:
    grid = twopanel(fig)
    clrmap    = get_color_map("viridis", over='r', under='k')
    
    maxerr = max((np.max(d[8]) for x0,x1,d in dats)) * 1.001
    for (x0,x1,dat),ax in zip(dats, grid):
      err = np.maximum(dat[8], 1e-20) 
      im = plot_color(np.log10(err), x0, x1, 
                      vmin=-16, vmax=log10(maxerr),  
                      cmap=clrmap, bar=False, axes=ax)
      rnd = dat[10]
      z1d = dat[0][:,0]
      b1d = dat[1][0,:]
      lvl = 1.0
      if np.max(rnd) > lvl:
        ax.contour(np.log10(z1d), np.log10(b1d), rnd.T, 
                   levels=[lvl], colors=['k'])
      ax.set_aspect('auto')
      ax.set_ylabel(r'$\log_{10}(B / \sqrt{D})$')
    grid[1].set_xlabel(r'$\log_{10}(z[c])$')
    color_bar(im, cax=ax.cax, label=r'$\log_{10}$(relative error)', extend='both')
  #
#

def load_data(path, eos):
  r = lambda t : load_grid(os.path.join(path, t % eos), 11)
  files = ["acc_eos%s_z_b_cold.dat", "acc_eos%s_z_b_hot.dat",
           "acc_eos%s_z_eps_Bzero.dat", "acc_eos%s_z_eps_Blarge.dat"]
  return map(r, files)
#


def main():
  path = sys.argv[1]
  
  for eos in ['ig', 'hyb']:
    dat_c, dat_h, dat_bz, dat_bl = load_data(path, eos)
    plot_acc_zeps(path, eos, dat_bz, dat_bl)
    plot_acc_zb(path, eos, dat_c, dat_h)
  #
#

main()
