#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
saturated_star.py
-----------------

Need to find a better target. Eclipse is too deep,
so it's hard to see the CDPP improvement for this one.

'''

import everest
import k2plr
import matplotlib.pyplot as pl
import numpy as np

star = everest.Everest(202063160)
inds = np.where(star.time >= 1939.11)[0]
time = star.time[inds]
flux = star.flux[inds]
fraw = star.fraw[inds]

v1_time, v1_flux = star.get_pipeline('everest1')
sff_time, sff_flux = star.get_pipeline('k2sff')
sff_flux *= np.nanmedian(flux)

fig, ax = pl.subplots(4, sharex = True, sharey = True)
ax[0].plot(time, fraw, 'k.', alpha = 0.3, ms = 3)
ax[1].plot(v1_time, v1_flux, 'k.', alpha = 0.3, ms = 3)
ax[2].plot(sff_time, sff_flux, 'k.', alpha = 0.3, ms = 3)
ax[3].plot(time, flux, 'k.', alpha = 0.3, ms = 3)

ax[0].margins(0.01, None)
pl.show()