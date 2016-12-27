#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
ridge_reg.py
------------

'''

from everest.config import EVEREST_DAT, EVEREST_SRC
from everest.missions.k2.aux import GetK2Campaign
from everest.utils import sort_like
import os
import matplotlib.pyplot as pl
from matplotlib.lines import Line2D
import numpy as np

campaign = 6.0
model = 'sPLD'

# Get statistics
outfile = os.path.join(EVEREST_SRC, 'missions', 'k2', 'tables', 'c%02d_%s.cdpp' % (int(campaign), model))
epic, kp, cdpp6r, cdpp6, cdpp6v, _, out, tot, saturated = np.loadtxt(outfile, unpack = True, skiprows = 2)
epic = np.array(epic, dtype = int)
saturated = np.array(saturated, dtype = int)

# Get only stars in this subcampaign
sub = np.array([s[0] for s in GetK2Campaign(campaign)], dtype = int)
inds = np.array([e in sub for e in epic])
epic = epic[inds]
kp = kp[inds]
cdpp6 = cdpp6[inds]
saturated = saturated[inds]
unsat = np.where(saturated == 0)

# Get only unsaturated stars
epic = epic[unsat]
kp = kp[unsat]
cdpp6 = cdpp6[unsat]

# Get the comparison model stats
epic_1, cdpp6_1 = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2', 
                             'tables', 'c%02d_everest1.cdpp' % int(campaign)), unpack = True)
cdpp6_1 = sort_like(cdpp6_1, epic, epic_1)  

# Plot
fig, axes = pl.subplots(1, 2, figsize = (10, 4))
alpha = min(0.1, 2000. / (1 + len(epic)))
y = (cdpp6 - cdpp6_1) / cdpp6_1
for ax in axes:
  ax.scatter(kp, y, color = 'b', marker = '.', alpha = alpha, zorder = -1)
  ax.set_xlim(11,18)
  ax.axhline(0, color = 'gray', lw = 2, zorder = -99, alpha = 0.5)
  ax.axhline(0.5, color = 'gray', ls = '--', lw = 2, zorder = -99, alpha = 0.5)
  ax.axhline(-0.5, color = 'gray', ls = '--', lw = 2, zorder = -99, alpha = 0.5)
  bins = np.arange(10.5,18.5,0.5)
  # Bin the CDPP
  by = np.zeros_like(bins) * np.nan
  for b, bin in enumerate(bins):
    i = np.where((y > -np.inf) & (y < np.inf) & (kp >= bin - 0.5) & (kp < bin + 0.5))[0]
    if len(i) > 10:
      by[b] = np.median(y[i])
  ax.plot(bins, by, 'k-', lw = 2)
  ax.set_xlabel('Kepler Magnitude', fontsize = 18)
  ax.set_rasterization_zorder(0)
axes[0].set_ylabel(r'$\frac{\mathrm{CDPP}_{\mathrm{L2}} - \mathrm{CDPP}_{\mathrm{EVEREST1}}}{\mathrm{CDPP}_{\mathrm{EVEREST1}}}$', fontsize = 22)
axes[0].set_ylim(-1,1)
axes[1].set_ylim(-0.3,0.3)

# Mark inset
axes[1].yaxis.tick_right()
axes[0].axhline(0.3, color = 'k', ls = '-', alpha = 1, lw = 1)
axes[0].axhline(-0.3, color = 'k', ls = '-', alpha = 1, lw = 1)
transFigure = fig.transFigure.inverted()
coord1 = transFigure.transform(axes[0].transData.transform([18,0.3]))
coord2 = transFigure.transform(axes[1].transData.transform([11,0.3]))
line1 = Line2D((coord1[0],coord2[0]),(coord1[1],coord2[1]),
                                transform=fig.transFigure, color = 'k', ls = '-', alpha = 1, lw = 1)
coord1 = transFigure.transform(axes[0].transData.transform([18,-0.3]))
coord2 = transFigure.transform(axes[1].transData.transform([11,-0.3]))
line2 = Line2D((coord1[0],coord2[0]),(coord1[1],coord2[1]),
                                transform=fig.transFigure, color = 'k', ls = '-', alpha = 1, lw = 1)
fig.lines = line1, line2

# Save
fig.savefig('ridge_reg.pdf', bbox_inches = 'tight')
pl.close()