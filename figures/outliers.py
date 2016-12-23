#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
outliers.py
-----------

'''

from everest.config import EVEREST_DAT, EVEREST_SRC
from everest.missions.k2.aux import GetK2Campaign
from everest.utils import sort_like
import os
import matplotlib.pyplot as pl
import numpy as np

fig, axes = pl.subplots(3, 3, figsize = (15, 10))
fig.subplots_adjust(left = 0.05, right = 0.975, bottom = 0.05, 
                    top = 0.95, hspace = 0.3, wspace = 0.20)
axes = axes.flatten()

# Loop over all campaigns
for campaign, ax in enumerate(axes):
  
  # DEBUG
  if campaign == 5 or campaign == 6:
    continue
  
  # Everest 2
  _,_,_,_,_,_, out_everest2, tot_everest2, _ = np.loadtxt(os.path.join(EVEREST_SRC, 'missions',
                                               'k2', 'tables', 'c%02d_nPLD.cdpp' % int(campaign)), 
                                               unpack = True, skiprows = 2)
  out_everest2 = np.array(out_everest2, dtype = int)
  tot_everest2 = np.array(tot_everest2, dtype = int)
  good_everest2 = tot_everest2 - out_everest2

  # Everest 1
  _, out_everest1, tot_everest1 = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2', 
                                  'tables', 'c%02d_everest1.out' % int(campaign)), 
                                  unpack = True, dtype = int)
  good_everest1 = tot_everest1 - out_everest1

  # K2SFF
  _, out_k2sff, tot_k2sff = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2', 
                            'tables', 'c%02d_k2sff.out' % int(campaign)), 
                            unpack = True, dtype = int)
  good_k2sff = tot_k2sff - out_k2sff

  # K2SC
  try:
    _, out_k2sc, tot_k2sc = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2', 
                            'tables', 'c%02d_k2sc.out' % int(campaign)), 
                            unpack = True, dtype = int)
    good_k2sc = tot_k2sc - out_k2sc
  except:
    good_k2sc = []
    
  # Get hist bounds
  all = np.concatenate([good_everest2, good_everest1, good_k2sff, good_k2sc])
  all = all[all > 0]
  i = np.argsort(all)
  a = int(0.01 * len(all))
  b = int(0.99 * len(all))
  xmin, xmax = all[i][a], all[i][b]
  xpad = 0.1 * (xmax - xmin)
  xmin -= xpad
  xmax += xpad
  
  # Plot
  ax.hist(good_k2sff, 25, range = (xmin, xmax), histtype = 'step', color = '#555555', label = 'K2SFF')
  ax.hist(good_k2sff, 25, range = (xmin, xmax), histtype = 'stepfilled', color = '#555555', alpha = 0.05, ec = None)
  if len(good_k2sc):
    ax.hist(good_k2sc, 25, range = (xmin, xmax), histtype = 'step', color = '#ee7600', label = 'K2SC')
    ax.hist(good_k2sc, 25, range = (xmin, xmax), histtype = 'stepfilled', color = '#ee7600', alpha = 0.05, ec = None)
  ax.hist(good_everest1, 25, range = (xmin, xmax), histtype = 'step', color = 'r', label = 'Everest1')
  ax.hist(good_everest1, 25, range = (xmin, xmax), histtype = 'stepfilled', color = 'r', alpha = 0.05, ec = None)
  ax.hist(good_everest2, 25, range = (xmin, xmax), histtype = 'step', color = '#0033ee', label = 'Everest2')
  ax.hist(good_everest2, 25, range = (xmin, xmax), histtype = 'stepfilled', color = '#0033ee', alpha = 0.05, ec = None)
  ax.margins(0, None)
  ax.legend(loc = 'upper left', fontsize = 10)
  ymax = ax.get_ylim()[1]
  ax.set_ylim(0, 1.2 * ymax)
  ax.set_title('C%02d' % campaign, fontsize = 16)
  
pl.show()