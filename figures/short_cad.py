#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
short_cad.py
------------

'''

import everest
import matplotlib.pyplot as pl
import numpy as np
from scipy.signal import savgol_filter

# Get the SC data
sc = everest.Everest(201601162, cadence = 'sc')

# Get the LC data
lc = everest.Everest(201601162, cadence = 'lc')

fig, ax = pl.subplots(3, figsize = (14, 8), sharex = True, sharey = True)
ax[0].plot(sc.time, sc.fraw, 'r.', alpha = 0.1, ms = 2)
ax[1].plot(sc.time, sc.flux, 'k.', alpha = 0.1, ms = 2)
ax[2].plot(sc.time, sc.flux, 'k.', alpha = 0.1, ms = 2)
ax[2].plot(lc.time, lc.flux, 'r.', alpha = 0.3, ms = 3)
pl.show()