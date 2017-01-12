#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
crowding.py
-----------

'''

import everest
import matplotlib.pyplot as pl
import numpy as np
import george
import k2plr
from scipy.ndimage import zoom

# For plotting the aperture
def PadWithZeros(vector, pad_width, iaxis, kwargs):
  vector[:pad_width[0]] = 0
  vector[-pad_width[1]:] = 0
  return vector

# Our targets
epics =   [202072978, 202103762, 218803648, 202733088, 211685048]
periods = [4.8840812, 1.3264960,   14.0532, 3.4338287, 0.7690710]
t0s =     [0.7787,         0.15,     -0.84,     0.675, 2306.8526]
durs =    [0.9,            0.15,       0.8,      0.35,       0.2]

# Set up plot
fig, ax = pl.subplots(5, 3, figsize = (7, 15))
for axis in ax.flatten():
  axis.set_xticklabels([])
  axis.set_yticklabels([])
ax[0,0].set_title('Everest 1.0', fontsize = 16)
ax[0,1].set_title('Everest 2.0', fontsize = 16)
plasma = pl.get_cmap('gray_r')
plasma.set_bad(alpha = 0)
  
for n, epic, period, t0, dur in zip(range(len(epics)), epics, periods, t0s, durs):
  
  # Get the data
  e2 = everest.Everest(epic)
  e1 = k2plr.EVEREST(epic, version = 1)

  # Whiten with the GP
  _, amp, tau = e2.kernel_params
  # HACK: GP for these stars is too strong, washes out eclipse
  if epic in [202072978, 218803648, 202733088]:
    amp /= 1000
  gp = george.GP(amp ** 2 * george.kernels.Matern32Kernel(tau ** 2))
  
  # Everest 2
  mask = []
  t0 += np.ceil((e2.time[0] - dur - t0) / period) * period
  for t in np.arange(t0, e2.time[-1] + dur, period):
    mask.extend(np.where(np.abs(e2.time - t) < dur / 2.)[0])
  e2_mask = np.array(list(set(np.concatenate([np.where(np.isnan(e2.flux))[0], mask]))))
  gp.compute(np.delete(e2.time, e2_mask), np.ones_like(np.delete(e2.time, e2_mask)) * np.nanmedian(e2.fraw_err))
  med = np.nanmedian(np.delete(e2.flux, e2_mask))
  y, _ = gp.predict(np.delete(e2.flux, e2_mask) - med, e2.time)
  fwhite2 = (e2.flux - y)
  fwhite2 /= np.nanmedian(fwhite2)
  tfold2 = (e2.time - t0 - period / 2.) % period - period / 2. 
  
  # Everest 1
  mask = []
  t0 += np.ceil((e1.time[0] - dur - t0) / period) * period
  for t in np.arange(t0, e1.time[-1] + dur, period):
    mask.extend(np.where(np.abs(e1.time - t) < dur / 2.)[0])
  e1.mask = np.array(list(set(np.concatenate([np.where(np.isnan(e1.flux))[0], mask]))))
  gp.compute(np.delete(e1.time, e1.mask), np.ones_like(np.delete(e1.time, e1.mask)) * np.nanmedian(e2.fraw_err))
  med = np.nanmedian(np.delete(e1.flux, e1.mask))
  y, _ = gp.predict(np.delete(e1.flux, e1.mask) - med, e1.time)
  fwhite1 = (e1.flux - y)
  fwhite1 /= np.nanmedian(fwhite1)
  tfold1 = (e1.time - t0 - period / 2.) % period - period / 2.
    
  # HACK: C0 data before t = 1937 has a small offset in the TPF timestamps (?!), 
  # leading to misfolded eclipses/transits. Let's exclude it.
  inds = np.where((np.abs(tfold2) < dur) & (e2.time > 1937.))[0]
  x2 = tfold2[inds]
  y2 = fwhite2[inds]
  inds = np.where((np.abs(tfold1) < dur) & (e1.time > 1937.))[0]
  x1 = tfold1[inds]
  y1 = fwhite1[inds]

  # Plot
  ax[n,0].plot(x1, y1, 'k.', alpha = 0.5, ms = 2)
  ax[n,1].plot(x2, y2, 'k.', alpha = 0.5, ms = 2)
  ax[n,0].set_ylabel('%d' % epic, fontsize = 14)
  
  # Get ylims
  yfin = np.delete(y2, np.where(np.isnan(y2)))
  lo, hi = yfin[np.argsort(yfin)][[25,-25]]
  pad = (hi - lo) * 0.2
  ylim = (lo - pad, hi + pad)
  ax[n,0].set_ylim(*ylim)
  ax[n,1].set_ylim(*ylim)
  ax[n,1].set_xlim(*ax[n,0].get_xlim())
  
  # Plot Hi Res image
  ny, nx = e2.pixel_images[0].shape
  contour = np.zeros((ny,nx))
  contour[np.where(e2.aperture)] = 1
  contour = np.lib.pad(contour, 1, PadWithZeros)
  highres = zoom(contour, 100, order = 0, mode='nearest') 
  extent = np.array([-1, nx, -1, ny])
  ax[n,2].imshow(e2.hires, aspect = 'auto', extent = (-0.5, nx - 0.5, -0.5, ny - 0.5), interpolation = 'bicubic', cmap = plasma)
  ax[n,2].contour(highres, levels=[0.5], extent=extent, origin='lower', colors='r', linewidths=1)
  ax[n,2].set_xticks([])
  ax[n,2].set_yticks([])
  ax[n,2].set_xlim(-0.7, nx - 0.3)
  ax[n,2].set_ylim(-0.7, ny - 0.3)
  for source in e2.nearby:
    ax[n,2].annotate('%.1f' % source['mag'], 
                     xy = (source['x'] - source['x0'], source['y'] - source['y0']), 
                     ha = 'center', va = 'center', size = 14, color = 'r', fontweight = 'bold') 

# I had to manually crop the apertures, so I plotted the
# figure interactively and saved it from there
pl.show()