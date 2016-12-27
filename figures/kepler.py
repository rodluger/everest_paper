#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
kepler.py
---------

'''

from everest.config import EVEREST_DAT, EVEREST_SRC
from everest.missions.k2.aux import GetK2Campaign
from everest.utils import sort_like
import os
import matplotlib.pyplot as pl
from matplotlib.lines import Line2D
import numpy as np

# Set up
fig, axes = pl.subplots(3, 3, figsize = (13, 8))
fig.subplots_adjust(left = 0.1, right = 0.975, bottom = 0.1, 
                    top = 0.95, hspace = 0.3, wspace = 0.20)
axes = axes.flatten()
axes[7].set_xlabel('Kepler Magnitude', fontsize = 24, labelpad = 20)
axes[3].set_ylabel(r'CDPP (ppm)', fontsize = 26, labelpad = 20)

# Get Kepler stats
kic, kepler_kp, kepler_cdpp6 = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2', 
                                         'tables', 'kepler.cdpp'), unpack = True) 

# Loop over all campaigns
for campaign, ax in enumerate(axes):
  
  # Get statistics
  outfile = os.path.join(EVEREST_SRC, 'missions', 'k2', 'tables', 'c%02d_nPLD.cdpp' % (int(campaign)))
  epic, kp, cdpp6r, cdpp6, cdpp6v, _, out, tot, saturated = np.loadtxt(outfile, unpack = True, skiprows = 2)
  epic = np.array(epic, dtype = int)
  saturated = np.array(saturated, dtype = int)
  unsat = np.where(saturated == 0)
  sat = np.where(saturated == 1)
  alpha = min(0.1, 500. / (1 + len(unsat[0])))
  
  # HACK: Campaign 0 magnitudes are reported only to the nearest tenth,
  # so let's add a little noise to spread them out for nicer plotting
  if campaign == 0:
    kp = kp + 0.1 * (0.5 - np.random.random(len(kp)))
      
  bins = np.arange(7.5,18.5,0.5) 
  ax.scatter(kepler_kp, kepler_cdpp6, color = 'y', marker = '.', alpha = 0.03, zorder = -1)
  ax.scatter(kp, cdpp6, color = 'b', marker = '.', alpha = alpha, zorder = -1)
  for x, y, color in zip([kepler_kp, kp], [kepler_cdpp6, cdpp6], ['#ffff66', '#6666ff']):
    by = np.zeros_like(bins) * np.nan
    for b, bin in enumerate(bins):
      i = np.where((y > -np.inf) & (y < np.inf) & (x >= bin - 0.5) & (x < bin + 0.5))[0]
      if len(i) > 10:
        by[b] = np.median(y[i])
    ax.plot(bins, by, 'o', color = color)
  ax.set_ylim(-10, 500)
  ax.set_xlim(8, 18)
  ax.set_title('C%02d' % campaign, fontsize = 18)
  ax.set_rasterization_zorder(0)

axes[0].scatter([-10,-10],[-10,-10],color='y',label='Kepler')
axes[0].scatter([-10,-10],[-10,-10],color='b',label='Everest')
axes[0].legend(loc = 'upper left', fontsize = 14)  
fig.savefig('cdpp_kepler.pdf', bbox_inches = 'tight')
pl.close()