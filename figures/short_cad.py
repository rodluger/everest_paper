#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
short_cad.py
------------

'''

import everest
from everest.math import Downbin, Interpolate
import matplotlib.pyplot as pl
import numpy as np

# Get the LC data
lc = everest.Everest(201601162, cadence = 'lc')
offset = 100 * (lc.time - 2017.69)

# Get the SC data
sc = everest.Everest(201601162, cadence = 'sc')
scflux = Interpolate(sc.time, sc.mask, sc.flux)
scflux_downbin = Downbin(scflux, len(lc.flux))

# Plot
fig, ax = pl.subplots(3, figsize = (10, 9))
ax[0].plot(sc.time, sc.fraw, 'k.', alpha = 0.1, ms = 2, zorder = -1)
ax[0].set_rasterization_zorder(0)
ax[1].plot(sc.time, sc.flux, 'k.', alpha = 0.1, ms = 2, zorder = -1)
ax[1].set_rasterization_zorder(0)
ax[2].plot(lc.time, scflux_downbin, 'k.', alpha = 0.3, ms = 3)
ax[2].plot(lc.time, lc.flux + offset, 'r.', alpha = 0.3, ms = 3)

# Limits
for axis in ax:
  axis.set_xlim(2017.69, 2057.34)
  axis.set_ylim(2.855e6, 2.951e6)

# Appearance
ax[2].set_xlabel('Time (BJD - 2454833)', fontsize = 28)
ax[0].set_ylabel('Raw', fontsize = 20)
ax[0].annotate('106.96 ppm', xy = (0.015, 0.95), color = 'k', fontsize = 16, 
               xycoords = 'axes fraction', ha = 'left', va = 'top')
ax[1].set_ylabel('EVEREST', fontsize = 20)
ax[1].annotate('13.85 ppm', xy = (0.015, 0.95), color = 'k', fontsize = 16, 
               xycoords = 'axes fraction', ha = 'left', va = 'top')
ax[2].set_ylabel('Downbinned', fontsize = 20)
ax[2].annotate('13.85 ppm', xy = (2021.23, 2.91668e6), color = 'k', fontsize = 16, xytext = (25, 30),
               textcoords = 'offset points', xycoords = 'data', ha = 'center',
               va = 'center', arrowprops = dict(arrowstyle = "-|>", color = 'k'))
ax[2].annotate('18.41 ppm', xy = (2021.23, 2.91161e6), color = 'r', fontsize = 16, xytext = (-10, -35),
               textcoords = 'offset points', xycoords = 'data', ha = 'center',
               va = 'center', arrowprops = dict(arrowstyle = "-|>", color = 'r'))
for tick in ax[0].get_xticklabels() + ax[1].get_xticklabels() + ax[2].get_xticklabels():
  tick.set_fontsize(14)
  
# Legend
ax[2].plot([0,1],[0,1],'k.',label='Short cadence')
ax[2].plot([0,1],[0,1],'r.',label='Long cadence')
ax[2].legend(loc = 'lower right', numpoints = 4, fontsize = 14)

# Save
fig.savefig('short_cad.pdf', bbox_inches = 'tight')