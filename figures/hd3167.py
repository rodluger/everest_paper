import everest
from scipy.signal import savgol_filter
from everest.mathutils import Interpolate, Smooth
from everest.missions.k2 import CDPP
import matplotlib.pyplot as pl
import numpy as np

# Planet params
EPIC = 220383386

# Set up the figure
fig, ax = pl.subplots(2, sharex=True, sharey=True, figsize=(13, 9))
fig.subplots_adjust(hspace=0.05)

# Everest
star = everest.Everest(EPIC)
time = star.time
med = np.nanmedian(star.flux)
baseline = Interpolate(star.time, star.mask, star.flux)
baseline = savgol_filter(baseline, 49, 2)
flux = star.flux - baseline + med
ax[1].plot(time, flux, 'k.', alpha=0.75, ms=3)

# K2SFF
time_sff, flux_sff = everest.k2.pipelines.get(EPIC, 'k2sff')
med = np.nanmedian(flux_sff)
ys = flux_sff - Smooth(flux_sff, 50)
M = np.nanmedian(ys)
MAD = 1.4826 * np.nanmedian(np.abs(ys - M))
out = []
for i, _ in enumerate(flux_sff):
    if (ys[i] > M + 3 * MAD) or (ys[i] < M - 3 * MAD):
        out.append(i)
out = np.array(out, dtype=int)
baseline = Interpolate(time_sff, out, flux_sff)
baseline = savgol_filter(baseline, 49, 2)
flux_sff = flux_sff - baseline + med
ax[0].plot(time_sff, flux_sff, 'k.', alpha=0.75, ms=3)

# CDPP
cdppe = star.cdpp
cdppf = CDPP(flux_sff)

# Appearance
ax[1].set_ylabel('EVEREST', fontsize=26)
ax[0].set_ylabel('K2SFF', fontsize=26)
ax[1].set_xlabel('Time (BJD - 2454833)', fontsize=28)
ax[1].annotate("%.2f ppm" % cdppe, xy=(0.025, 0.96), xycoords="axes fraction",
               ha="left", va="top", fontsize=18, bbox=dict(fc="w", ec="1"),
               color="r")
ax[0].annotate("%.2f ppm" % cdppf, xy=(0.025, 0.96), xycoords="axes fraction",
               ha="left", va="top", fontsize=18, bbox=dict(fc="w", ec="1"),
               color="r")
ax[0].margins(0, None)
ax[0].set_ylim(0.99875 * med, 1.00125 * med)
for tick in ax[0].get_xticklabels() + ax[1].get_xticklabels():
    tick.set_fontsize(16)

# Save!
fig.savefig("hd3167.pdf", bbox_inches="tight")
