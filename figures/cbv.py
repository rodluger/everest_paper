#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
cbv.py
------

'''

from everest.missions.k2.sysrem import GetCBVs, GetChunk
from everest.config import EVEREST_DAT
import os
import matplotlib.pyplot as pl
import numpy as np

# Let's do C02 as an example
campaign = 2
clobber = False

# We're going to plot the CBVs on the CCD
fig = [None, None, None]
ax = [None, None, None]
for n in range(1, 3):
  fig[n], ax[n] = pl.subplots(5, 5, figsize = (9, 9))
  fig[n].subplots_adjust(wspace = 0.025, hspace = 0.025)
  ax[n] = [None] + list(ax[n].flatten())
  for axis in [ax[n][1], ax[n][5], ax[n][21], ax[n][25]]:
    axis.set_visible(False)
  for i in range(1, 25):
    ax[n][i].set_xticks([])
    ax[n][i].set_yticks([])
    ax[n][i].annotate('%02d' % i, (0.5, 0.5), 
                      va = 'center', ha = 'center',
                      xycoords = 'axes fraction',
                      color = 'k', fontsize = 60, alpha = 0.05)
    ax[n][i].margins(0.1, 0.1)
  
# Get the CBVs
for module in range(2, 25):
  X = GetCBVs(campaign, module = module, model = 'nPLD', clobber = clobber)
  if X is not None:
  
    # Get the timeseries info
    infofile = os.path.join(EVEREST_DAT, 'k2', 'cbv', 'c%02d' % campaign, str(module), 'nPLD', 'info.npz')
    info = np.load(infofile)
    time = info['time']
    nstars = info['nstars']
    breakpoints = info['breakpoints']
    
    # Plot the CBVs
    for b, color in zip(range(len(breakpoints)), ['b', 'r']):
      inds = GetChunk(time, breakpoints, b)
      for n in range(1, 3):
        if n == 2 and module == 20:
          alpha = 0.25
        else:
          alpha = 1.0
        ax[n][module].plot(time[inds], X[inds,n], color = color, lw = 2, alpha = alpha)
          
# Save the figures
for n in range(1, 3):
  fig[n].savefig('cbv%d.pdf' % n, bbox_inches = 'tight')
  pl.close(fig[n])