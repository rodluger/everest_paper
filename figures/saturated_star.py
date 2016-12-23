#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
saturated_star.py
-----------------


'''

import everest
import k2plr
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import numpy as np
from scipy.ndimage import zoom

# Everest 2.0
star = everest.Everest(202063160)
inds = np.where(star.time >= 1959)[0]
time = star.time[inds]
flux = star.flux[inds]
flux /= np.nanmedian(flux)

# Raw
fraw = star.fraw[inds]
fraw /= np.nanmedian(fraw)

# Everest 1.0
v1_time, v1_flux = star.get_pipeline('everest1')
inds = np.where(v1_time >= 1959)[0]
v1_time = v1_time[inds]
v1_flux = v1_flux[inds]
v1_flux /= np.nanmedian(v1_flux)

# Plot
fig = pl.figure(figsize = (10, 9))
ax_raw = pl.subplot2grid((120, 120), (0,  0), colspan=80, rowspan=35)
ax_v1 = pl.subplot2grid((120, 120), (40,  0), colspan=80, rowspan=35, sharex = ax_raw, sharey = ax_raw)
ax_v2 = pl.subplot2grid((120, 120), (80,  0), colspan=80, rowspan=35, sharex = ax_raw, sharey = ax_raw)
ax_ap = pl.subplot2grid((120, 120), (0,  85), colspan=35, rowspan=115)

# Light curves
ax_raw.plot(time, fraw, 'k.', alpha = 0.7, ms = 3)
ax_raw.annotate('%.2f ppm' % star.cdppr, xy = (0.025, 0.95), xycoords = 'axes fraction',
                ha = 'left', va = 'top', color = 'r', fontsize = 12)
ax_v1.plot(v1_time, v1_flux, 'k.', alpha = 0.7, ms = 3)
ax_v1.annotate('%.2f ppm' % 9.38, xy = (0.025, 0.95), xycoords = 'axes fraction',
                ha = 'left', va = 'top', color = 'r', fontsize = 12)
ax_v2.plot(time, flux, 'k.', alpha = 0.7, ms = 3)
ax_v2.annotate('%.2f ppm' % star.cdpp, xy = (0.025, 0.95), xycoords = 'axes fraction',
                ha = 'left', va = 'top', color = 'r', fontsize = 12)
ax_raw.margins(0.01, 0.1)
ax_raw.get_yaxis().set_major_locator(MaxNLocator(4))
ax_v1.get_yaxis().set_major_locator(MaxNLocator(4))
ax_v2.get_yaxis().set_major_locator(MaxNLocator(4))
ax_v2.set_xlabel('Time (BJD - 2454833)', fontsize = 18)
#ax_raw.set_xticklabels([])
#ax_v1.set_xticklabels([])
ax_raw.set_ylabel('Raw', fontsize = 18)
ax_v1.set_ylabel('Everest 1', fontsize = 18)
ax_v2.set_ylabel('Everest 2', fontsize = 18)

# Aperture
plasma = pl.get_cmap('Greys')
plasma.set_bad(alpha = 0)
def PadWithZeros(vector, pad_width, iaxis, kwargs):
  vector[:pad_width[0]] = 0
  vector[-pad_width[1]:] = 0
  return vector
ny, nx = star.pixel_images[2].shape
contour = np.zeros((ny,nx))
contour[np.where(star.aperture)] = 1
contour = np.lib.pad(contour, 1, PadWithZeros)
highres = zoom(contour, 100, order = 0, mode='nearest') 
extent = np.array([-1, nx, -1, ny])
ax_ap.imshow(star.pixel_images[2], aspect = 'auto', interpolation = 'nearest', cmap = plasma)
ax_ap.contour(highres, levels=[0.5], extent=extent, origin='lower', colors='r', linewidths=1)
ax_ap.axvspan(14.5, 17.5, color = 'r', alpha = 0.1, ec = 'none')
ax_ap.set_xticks([])
ax_ap.set_yticks([])       
ax_ap.set_xlim(6, 24)
ax_ap.set_ylim(10, 28)

# Save
fig.savefig('saturated_star.pdf', bbox_inches = 'tight')