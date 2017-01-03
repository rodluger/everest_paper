#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
201345483.py
------------

'''

import everest
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import numpy as np
import george

# Planet params
EPIC = 201345483
t0 = 1978.25
period = 1.72925
dur = 0.1

# Get the Everest 2 data
evr2 = everest.Everest(EPIC)
evr2.mask_planet(t0, period, dur)
evr2.compute()

# Get the K2SFF data
k2sff_time, k2sff_flux = everest.k2.pipelines.get(EPIC, 'k2sff')

# Fill in missing cadences
tol = 0.005
if not ((len(k2sff_time) == len(evr2.time)) and (np.abs(k2sff_time[0] - evr2.time[0]) < tol) and (np.abs(k2sff_time[-1] - evr2.time[-1]) < tol)):
  ftmp = np.zeros_like(evr2.time)
  j = 0
  for i, t in enumerate(evr2.time):
    if np.abs(k2sff_time[j] - t) < tol:
      ftmp[i] = k2sff_flux[j]
      j += 1
      if j == len(k2sff_time) - 1:
        break
  k2sff_flux = ftmp
  k2sff_time = np.array(evr2.time)

# Set up the plot
fig = pl.figure(figsize = (13, 9))
axk2sff = pl.subplot2grid((2, 10), (0, 0), colspan = 7, rowspan = 1)
axk2sfff = pl.subplot2grid((2, 10), (0, 7), colspan = 3, rowspan = 1)
axevr2 = pl.subplot2grid((2, 10), (1, 0), colspan = 7, rowspan = 1)
axevr2f = pl.subplot2grid((2, 10), (1, 7), colspan = 3, rowspan = 1)

# Plot the K2SFF light curve
axk2sff.plot(k2sff_time, k2sff_flux, 'k.', alpha = 0.3, ms = 3)
axk2sff.set_ylim(9700, 10600)
axk2sff.set_xlim(1977.26, 2057.33)

# Plot the Everest 2 light curve
y = np.array(evr2.flux)
# [HACK] Remove the missing mid-campaign data
y[np.where((evr2.time >= 2014.98) & (evr2.time <= 2017.88))] = np.nan
axevr2.plot(evr2.time, y, 'k.', alpha = 0.3, ms = 3)
axevr2.set_ylim(9700, 10600)
axevr2.set_xlim(1977.26, 2057.33)

# Plot the whitened and folded K2SFF light curve
_, amp, tau = evr2.kernel_params
gp = george.GP(amp ** 2 * george.kernels.Matern32Kernel(tau ** 2))
gp.compute(evr2.apply_mask(k2sff_time), evr2.apply_mask(evr2.fraw_err))
med = np.nanmedian(evr2.apply_mask(k2sff_flux))
y, _ = gp.predict(evr2.apply_mask(k2sff_flux) - med, k2sff_time)
fwhite = (k2sff_flux - y)
fwhite /= np.nanmedian(fwhite)
tfold = (k2sff_time - t0 - period / 2.) % period - period / 2. 
inds = np.where(np.abs(tfold) < 1.5 * dur)[0]
x = tfold[inds]
y = fwhite[inds]
axk2sfff.plot(x, y, 'k.', alpha = 0.75, ms = 3)
axk2sfff.set_ylim(0.973, 1.005)
axk2sfff.yaxis.tick_right()
axk2sfff.get_yaxis().set_major_locator(MaxNLocator(4))
axk2sfff.get_xaxis().set_major_locator(MaxNLocator(4))

# Plot the whitened and folded Everest 2 light curve
_, amp, tau = evr2.kernel_params
gp = george.GP(amp ** 2 * george.kernels.Matern32Kernel(tau ** 2))
gp.compute(evr2.apply_mask(evr2.time), evr2.apply_mask(evr2.fraw_err))
med = np.nanmedian(evr2.apply_mask(evr2.flux))
y, _ = gp.predict(evr2.apply_mask(evr2.flux) - med, evr2.time)
fwhite = np.array(evr2.flux)
# [HACK] Remove the missing mid-campaign data
fwhite[np.where((evr2.time >= 2014.98) & (evr2.time <= 2017.88))] = np.nan
fwhite -= y
fwhite /= np.nanmedian(fwhite)
tfold = (evr2.time - t0 - period / 2.) % period - period / 2. 
inds = np.where(np.abs(tfold) < 1.5 * dur)[0]
x = tfold[inds]
y = fwhite[inds]
axevr2f.plot(x, y, 'k.', alpha = 0.75, ms = 3)
axevr2f.set_ylim(0.973, 1.005)
axevr2f.yaxis.tick_right()
axevr2f.get_yaxis().set_major_locator(MaxNLocator(4))
axevr2f.get_xaxis().set_major_locator(MaxNLocator(4))

# Labels and stuff
axevr2.set_xlabel('Time (BJD - 2454833)', fontsize = 28)
axevr2f.set_xlabel('Time (days)', fontsize = 28)
axevr2.set_ylabel('EVEREST', fontsize = 26)
axk2sff.set_ylabel('K2SFF', fontsize = 26)
axevr2.annotate('119.66 ppm', xy = (0.025, 0.96), xycoords = 'axes fraction', 
                ha = 'left', va = 'top', color = 'r', fontsize = 18)
axk2sff.annotate('285.45 ppm', xy = (0.025, 0.96), xycoords = 'axes fraction', 
                 ha = 'left', va = 'top', color = 'r', fontsize = 18)    
for tick in axevr2.get_xticklabels() + axevr2f.get_xticklabels() + axk2sff.get_xticklabels() + axk2sfff.get_xticklabels():
  tick.set_fontsize(16)
axevr2.get_yaxis().set_major_locator(MaxNLocator(6))
axk2sff.get_yaxis().set_major_locator(MaxNLocator(6)) 
            
# Save
fig.savefig('201345483.pdf', bbox_inches = 'tight')