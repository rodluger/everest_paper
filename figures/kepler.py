#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
kepler.py
---------

'''

from everest.config import EVEREST_DAT, EVEREST_SRC
from everest.missions.k2.utils import GetK2Campaign
from everest.utils import sort_like
import os
import matplotlib.pyplot as pl
from matplotlib.lines import Line2D
import numpy as np

# Set up
fig, axes = pl.subplots(4, 4, figsize=(13, 10))
fig.subplots_adjust(left=0.1, right=0.975, bottom=0.1,
                    top=0.95, hspace=0.3, wspace=0.20)
axes = axes.flatten()
axes[14].scatter([-10, -11, -12], [-10, -11, -12], color='y', label='Kepler')
axes[14].scatter([-10, -11, -12], [-10, -11, -12], color='b', label='Everest')
axes[14].set_xlim(0, 1)
axes[14].set_ylim(0, 1)
axes[14].legend(loc='center', fontsize=18, scatterpoints=3)
axes[14].axis('off')
axes[15].axis('off')
axes = axes[:-2]
fig.text(0.5, 0.04, 'Kepler Magnitude',
         fontsize=22, ha='center', va='center')
fig.text(0.02, 0.5, r'CDPP (ppm)', fontsize=24,
         ha='center', va='center', rotation='vertical')

# Get Kepler stats
kic, kepler_kp, kepler_cdpp6 = np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2',
                                                       'tables', 'kepler.cdpp'), unpack=True)

# Cumulative lists
kp_all = []
cdpp6_all = []

# Loop over all campaigns
campaigns = [0, 1, 2, 3, 4, 5, 6, 7, 8, 102, 111, 112, 12, 13]
for i, ax in enumerate(axes):

    campaign = campaigns[i]

    # Get statistics
    outfile = os.path.join(EVEREST_SRC, 'missions', 'k2',
                           'tables', 'c%02d_nPLD.cdpp' % (int(campaign)))
    epic, kp, cdpp6r, cdpp6, cdpp6v, _, out, tot, saturated = np.loadtxt(
        outfile, unpack=True, skiprows=2)
    epic = np.array(epic, dtype=int)
    saturated = np.array(saturated, dtype=int)
    unsat = np.where(saturated == 0)
    sat = np.where(saturated == 1)
    alpha = min(0.1, 500. / (1 + len(unsat[0])))

    # HACK: Campaign 0 magnitudes are reported only to the nearest tenth,
    # so let's add a little noise to spread them out for nicer plotting
    if campaign == 0:
        kp = kp + 0.1 * (0.5 - np.random.random(len(kp)))

    # Append to running lists
    kp_all.extend(list(kp))
    cdpp6_all.extend(list(cdpp6))

    bins = np.arange(7.5, 18.5, 0.5)
    ax.scatter(kepler_kp, kepler_cdpp6, color='y',
               marker='.', alpha=0.03, zorder=-1)
    ax.scatter(kp, cdpp6, color='b', marker='.', alpha=alpha, zorder=-1)
    for x, y, color in zip([kepler_kp, kp], [kepler_cdpp6, cdpp6], ['#ffff66', '#6666ff']):
        by = np.zeros_like(bins) * np.nan
        for b, bin in enumerate(bins):
            i = np.where((y > -np.inf) & (y < np.inf) &
                         (x >= bin - 0.5) & (x < bin + 0.5))[0]
            if len(i) > 10:
                by[b] = np.median(y[i])
        ax.plot(bins, by, 'o', color=color)
    ax.set_ylim(-10, 500)
    ax.set_xlim(8, 18)
    ax.set_title('C%02d' % campaign, fontsize=18)
    ax.set_rasterization_zorder(0)

fig.savefig('cdpp_kepler.pdf', bbox_inches='tight')
pl.close()

# Now plot the cumulative one
kp_all = np.array(kp_all)
cdpp6_all = np.array(cdpp6_all)
fig, ax = pl.subplots(1, 1, figsize=(7, 5))
ax.set_xlabel('Kepler Magnitude', fontsize=28, labelpad=20)
ax.set_ylabel(r'CDPP (ppm)', fontsize=28, labelpad=20)
bins = np.arange(7.5, 18.5, 0.5)
ax.scatter(kepler_kp, kepler_cdpp6, color='y',
           marker='.', alpha=0.03, zorder=-1)
ax.scatter(kp_all[::2], cdpp6_all[::2], color='b',
           marker='.', alpha=0.01, zorder=-1)
for x, y, color in zip([kepler_kp, kp_all], [kepler_cdpp6, cdpp6_all], ['#ffff66', '#6666ff']):
    by = np.zeros_like(bins) * np.nan
    for b, bin in enumerate(bins):
        i = np.where((y > -np.inf) & (y < np.inf) &
                     (x >= bin - 0.5) & (x < bin + 0.5))[0]
        if len(i) > 10:
            by[b] = np.median(y[i])
    ax.plot(bins, by, 'o', color=color)
ax.set_ylim(-10, 500)
ax.set_xlim(8, 18)
ax.set_rasterization_zorder(0)
ax.scatter([-10, -10], [-10, -10], color='y', label='Kepler')
ax.scatter([-10, -10], [-10, -10], color='b', label='Everest')
ax.legend(loc='upper left', fontsize=16, scatterpoints=3)
for tick in ax.get_xticklabels() + ax.get_yticklabels():
    tick.set_fontsize(16)
fig.savefig('cdpp_kepler_all.pdf', bbox_inches='tight')
pl.close()
