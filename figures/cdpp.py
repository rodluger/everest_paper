#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
cdpp.py
-------

'''

from everest.config import EVEREST_DAT, EVEREST_SRC
from everest.missions.k2.aux import GetK2Campaign
from everest.utils import sort_like
import os
import matplotlib.pyplot as pl
from matplotlib.lines import Line2D
import numpy as np

for model in ['everest1', 'k2sff', 'k2sc']:
  
  # Set up the campaign-by-campaign figure
  if model == 'k2sc':
    fig, axes = pl.subplots(2, 2, figsize = (9, 6))
    fig.subplots_adjust(left = 0.1, right = 0.975, bottom = 0.1, 
                        top = 0.95, hspace = 0.3, wspace = 0.20)
    axes = axes.flatten()
    axes = [None, None, None, axes[0], axes[1], axes[2], axes[3]]
    fig.text(0.5, 0.01, 'Kepler Magnitude', fontsize = 22, ha='center', va='center')
    fig.text(0.015, 0.5, r'$\frac{\mathrm{CDPP}_{\mathrm{EVEREST2}} - \mathrm{CDPP}_{\mathrm{%s}}}{\mathrm{CDPP}_{\mathrm{%s}}}$' % 
            (model.upper(), model.upper()), fontsize = 24, ha='center', va='center', rotation='vertical')

  else:
    fig, axes = pl.subplots(3, 3, figsize = (13, 10))
    fig.subplots_adjust(left = 0.1, right = 0.975, bottom = 0.1, 
                        top = 0.95, hspace = 0.3, wspace = 0.20)
    axes = axes.flatten()
    axes[7].set_xlabel('Kepler Magnitude', fontsize = 32, labelpad = 20)
    axes[3].set_ylabel(r'$\frac{\mathrm{CDPP}_{\mathrm{EVEREST2}} - \mathrm{CDPP}_{\mathrm{%s}}}{\mathrm{CDPP}_{\mathrm{%s}}}$' % 
                      (model.upper(), model.upper()), fontsize = 34, labelpad = 20)
  
  # Running lists
  yunsat = []
  kpunsat = []
  ysat = []
  kpsat = []
    
  # Loop over all campaigns
  for campaign, ax in enumerate(axes):
    
    # Some campaigns are missing from certain pipelines
    if ax is None:
      continue
    
    # Get statistics
    outfile = os.path.join(EVEREST_SRC, 'missions', 'k2', 'tables', 'c%02d_nPLD.cdpp' % (int(campaign)))
    epic, kp, cdpp6r, cdpp6, cdpp6v, _, out, tot, saturated = np.loadtxt(outfile, unpack = True, skiprows = 2)
    epic = np.array(epic, dtype = int)
    saturated = np.array(saturated, dtype = int)
    unsat = np.where(saturated == 0)
    sat = np.where(saturated == 1)
  
    # HACK: Campaign 0 magnitudes are reported only to the nearest tenth,
    # so let's add a little noise to spread them out for nicer plotting
    if campaign == 0:
      kp = kp + 0.1 * (0.5 - np.random.random(len(kp)))

    # Get the comparison model stats
    epic_1, cdpp6_1 = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2', 
                                 'tables', 'c%02d_%s.cdpp' % (int(campaign), model)), unpack = True)
    cdpp6_1 = sort_like(cdpp6_1, epic, epic_1)  
        
    # Plot
    alpha = min(0.1, 500. / (1 + len(epic)))
    y = (cdpp6 - cdpp6_1) / cdpp6_1
    ax.scatter(kp[unsat], y[unsat], color = 'b', marker = '.', alpha = alpha, zorder = -1)
    ax.scatter(kp[sat], y[sat], color = 'r', marker = '.', alpha = alpha, zorder = -1)
    ax.set_xlim(8,18)
    ax.axhline(0, color = 'gray', lw = 2, zorder = -99, alpha = 0.5)
    ax.axhline(0.5, color = 'gray', ls = '--', lw = 2, zorder = -99, alpha = 0.5)
    ax.axhline(-0.5, color = 'gray', ls = '--', lw = 2, zorder = -99, alpha = 0.5)
    bins = np.arange(7.5,18.5,0.5)

    # Bin the CDPP
    by = np.zeros_like(bins) * np.nan
    for b, bin in enumerate(bins):
      i = np.where((y > -np.inf) & (y < np.inf) & (kp >= bin - 0.5) & (kp < bin + 0.5))[0]
      if len(i) > 10:
        by[b] = np.median(y[i])
    ax.plot(bins[:9], by[:9], 'r--', lw = 2)
    ax.plot(bins[8:], by[8:], 'k-', lw = 2)
    ax.set_ylim(-1,1)
    ax.set_title('C%02d' % campaign, fontsize = 22)
    ax.set_rasterization_zorder(0)
    
    # Append to running lists
    yunsat.extend(y[unsat])
    ysat.extend(y[sat])
    kpunsat.extend(kp[unsat])
    kpsat.extend(kp[sat])
    
  fig.savefig('cdpp_%s.pdf' % model, bbox_inches = 'tight')
  pl.close()
  
  # Plot all campaigns at once
  fig, ax = pl.subplots(1, figsize = (9,6))
  ax.set_xlabel('Kepler Magnitude', fontsize = 22, labelpad = 20)
  ax.set_ylabel(r'$\frac{\mathrm{CDPP}_{\mathrm{EVEREST2}} - \mathrm{CDPP}_{\mathrm{%s}}}{\mathrm{CDPP}_{\mathrm{%s}}}$' % 
                (model.upper(), model.upper()), fontsize = 24, labelpad = 20)
  alpha = min(0.1, 1000. / (1 + len(yunsat)))
  
  ax.scatter(kpunsat, yunsat, color = 'b', marker = '.', alpha = alpha, zorder = -1)
  ax.scatter(kpsat, ysat, color = 'r', marker = '.', alpha = alpha, zorder = -1)
  ax.set_xlim(8,18)
  ax.axhline(0, color = 'gray', lw = 2, zorder = -99, alpha = 0.5)
  ax.axhline(0.5, color = 'gray', ls = '--', lw = 2, zorder = -99, alpha = 0.5)
  ax.axhline(-0.5, color = 'gray', ls = '--', lw = 2, zorder = -99, alpha = 0.5)
  bins = np.arange(7.5,18.5,0.5)

  # Bin the CDPP
  y = np.append(ysat, yunsat)
  kp = np.append(kpsat, kpunsat)
  by = np.zeros_like(bins) * np.nan
  for b, bin in enumerate(bins):
    i = np.where((y > -np.inf) & (y < np.inf) & (kp >= bin - 0.5) & (kp < bin + 0.5))[0]
    if len(i) > 10:
      by[b] = np.median(y[i])
  ax.plot(bins[:9], by[:9], 'r--', lw = 2)
  ax.plot(bins[8:], by[8:], 'k-', lw = 2)
  ax.set_ylim(-1,1)
  ax.set_rasterization_zorder(0)
  fig.savefig('cdpp_%s_all.pdf' % model, bbox_inches = 'tight')
  pl.close()
