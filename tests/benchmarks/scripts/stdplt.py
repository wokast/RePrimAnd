# -*- coding: utf-8 -*-
import os
import warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import ScalarFormatter
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import ImageGrid
from cycler import cycler
from math import *

def get_color_map(name, over=None, under=None, masked=None):
  cmap  = plt.get_cmap(name)
  if over:
    cmap.set_over(over)    
  if under:
    cmap.set_under(under)
  if masked:
    cmap.set_bad(masked)
  return cmap
#


def twopanel(fig):
  grid = ImageGrid(fig, 111,         
                 nrows_ncols=(2,1),
                 aspect=False,
                 axes_pad=0.1,
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="2%",
                 cbar_pad=0.1)
  return grid
#

def color_bar(image=None, cax=None, ax=None, label=None, ticks=None, **kwargs):
  if (image is None):
    cb  = plt.colorbar(cax=cax, ax=ax, **kwargs)
  else:      
    cb  = plt.colorbar(image, cax=cax, ax=ax, **kwargs)
  #
  if (label is not  None):
    cb.set_label(label)
  #
  if (ticks is not  None):
    cb.set_ticks(ticks)
    cb.set_ticklabels(ticks)
  #
  return cb
#


def plot_color(var, x0, x1, axes=None,  vmin=None, vmax=None, 
               cmap=None, interpolation='nearest',
               bar=False, barlabel=None, barticks=None, barshrink=1,
               barextend='neither',  **kwargs):
  x0, x1 = np.array(x0), np.array(x1)
  dx     = (x1-x0) / var.shape

  if (vmin is None):
    vmin=var.min()
  #
  if (vmax is None):
    vmax=var.max()
  #
  axes = plt.gca() if axes is None else axes
  
  ext   = [x0[0]-0.5*dx[0],
          x1[0]+0.5*dx[0],
          x0[1]-0.5*dx[1],
          x1[1]+0.5*dx[1]]
  z     = var.transpose()
  im    = axes.imshow(z, vmin=vmin, vmax=vmax, interpolation=interpolation, 
                   cmap=cmap, extent=ext, origin='lower', **kwargs)
  if bar:
    color_bar(im, ax=axes, shrink=barshrink, label=barlabel, ticks=barticks,
              extend=barextend)
  #
  
  return im
#


class Figure():
  def __init__(self, name, path='./', formats=['pdf'], 
               figwidth=1., figaspect=1.618, **kwargs):
    self.filename = os.path.join(path, name)
    self.formats = formats
    self.figargs = kwargs
    fw = 6.5*figwidth
    fh = fw/figaspect
    self.figargs.setdefault('figsize', [fw,fh])
  #
  def __enter__(self):
    self.fig = plt.figure(**self.figargs)
    return self.fig
  #
  def __exit__(self, exc, value, traceback):
    if exc is None:
      with warnings.catch_warnings():
        warnings.simplefilter("ignore") #annoying tight layout warning
        plt.tight_layout()
        for fm in self.formats:
          self.fig.savefig("%s.%s" % (self.filename, fm))
        #
      #
    else:
      print(value)
    #
    return True
  #
#



def load_grid(path, ncols):
  mkempty = lambda : [[] for i in range(ncols)]
  # [[]]*11 wont work here (-> 11 refs to same [])
  x1 = mkempty()
  x2 = mkempty()
  f = open(path, 'r')
  for l in f:
    c = l.split()
    if ((len(c) == 0) and (len(x1[0]) != 0)):
      for e1,e2 in zip(x1,x2):
        e2.append(e1)
      x1 = mkempty() 
      
    else:
      for e1,c1 in zip(x1,c):
        e1.append(float(c1))
      #
    #
  #
  x2 = [np.array(e) for e in x2]
  p0 = [log10(x2[0].min()), log10(x2[1].min())]
  p1 = [log10(x2[0].max()), log10(x2[1].max())]
  return p0, p1, x2
#
