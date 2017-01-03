#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
short_cad_stats.py
------------------

'''

from everest.config import EVEREST_DAT, EVEREST_SRC
from everest.missions.k2.aux import GetK2Campaign
from everest.utils import sort_like
import os
import matplotlib.pyplot as pl
import numpy as np

campaign = np.arange(9)
model = 'nPLD.sc'

# Running lists
xsat = []
ysat = []
xunsat = []
yunsat = []
xall = []
yall = []
epics = []

# Plot
for camp in campaign:

  # Load all stars
  sub = np.array(GetK2Campaign(camp, cadence = 'sc', epics_only = True), dtype = int)
  outfile = os.path.join(EVEREST_SRC, 'missions', 'k2', 'tables', 'c%02d_%s.cdpp' % (int(camp), model))
  epic, kp, cdpp6r, cdpp6, saturated = np.loadtxt(outfile, unpack = True, skiprows = 2)
  epic = np.array(epic, dtype = int)
  saturated = np.array(saturated, dtype = int)

  # Get only stars in this subcamp
  inds = np.array([e in sub for e in epic])
  epic = epic[inds]
  kp = kp[inds]
  # HACK: camp 0 magnitudes are reported only to the nearest tenth,
  # so let's add a little noise to spread them out for nicer plotting
  kp = kp + 0.1 * (0.5 - np.random.random(len(kp)))
  cdpp6r = cdpp6r[inds]
  cdpp6 = cdpp6[inds]
  saturated = saturated[inds]
  sat = np.where(saturated == 1)
  unsat = np.where(saturated == 0)
  if not np.any([not np.isnan(x) for x in cdpp6]):
    continue

  # Get the long cadence stats
  compfile = os.path.join(EVEREST_SRC, 'missions', 'k2', 'tables', 'c%02d_%s.cdpp' % (int(camp), model[:-3]))
  epic_1, _, _, cdpp6_1, _, _, _, _, saturated = np.loadtxt(compfile, unpack = True, skiprows = 2)
  epic_1 = np.array(epic_1, dtype = int)
  inds = np.array([e in sub for e in epic_1])
  epic_1 = epic_1[inds]
  cdpp6_1 = cdpp6_1[inds]
  cdpp6_1 = sort_like(cdpp6_1, epic, epic_1) 
  x = kp
  y = (cdpp6 - cdpp6_1) / cdpp6_1
  
  # Append to running lists
  xsat.extend(x[sat])
  ysat.extend(y[sat])
  xunsat.extend(x[unsat])
  yunsat.extend(y[unsat])
  xall.extend(x)
  yall.extend(y)
  epics.extend(epic)

# Plot the equivalent of the Aigrain+16 figure
fig, ax = pl.subplots(1, figsize = (9,6))
ax.set_xlabel('Kepler Magnitude', fontsize = 22, labelpad = 20)
ax.set_ylabel(r'$\frac{\mathrm{CDPP}_{\mathrm{SC}} - \mathrm{CDPP}_{\mathrm{LC}}}{\mathrm{CDPP}_{\mathrm{LC}}}$', 
              fontsize = 24, labelpad = 20)

# Scatter
ax.scatter(xunsat, yunsat, color = 'b', marker = '.', alpha = 0.35, zorder = -1, picker = True)
ax.scatter(xsat, ysat, color = 'r', marker = '.', alpha = 0.35, zorder = -1, picker = True)
ax.set_ylim(-1,1)
ax.set_xlim(8,18)
ax.axhline(0, color = 'gray', lw = 2, zorder = -99, alpha = 0.5)
ax.axhline(0.5, color = 'gray', ls = '--', lw = 2, zorder = -99, alpha = 0.5)
ax.axhline(-0.5, color = 'gray', ls = '--', lw = 2, zorder = -99, alpha = 0.5)

# Bin the CDPP
yall = np.array(yall)
xall = np.array(xall)
bins = np.arange(7.5,18.5,0.5)
by = np.zeros_like(bins) * np.nan
for b, bin in enumerate(bins):
  i = np.where((yall > -np.inf) & (yall < np.inf) & (xall >= bin - 0.5) & (xall < bin + 0.5))[0]
  if len(i) > 10:
    by[b] = np.median(yall[i])
ax.plot(bins[:9], by[:9], 'r--', lw = 2)
ax.plot(bins[8:], by[8:], 'k-', lw = 2)

for tick in ax.get_xticklabels() + ax.get_yticklabels():
  tick.set_fontsize(18)

# Save
fig.savefig('short_cad_stats.pdf', bbox_inches = 'tight')
pl.close()