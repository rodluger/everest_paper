#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
crossval.py
-----------

'''

import everest
from everest.config import EVEREST_DAT, EVEREST_SRC
from everest.utils import Formatter
from everest.math import Chunks
import os
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import numpy as np
import george

if not os.path.exists('crossval.npz'):

  # Get the model
  star = everest.nPLD(206103150, debug = True)

  # Let's just do the first chunk
  b = 0
  brkpt = star.breakpoints[b]

  # Let's just do the first PLD order
  star.lam = [[0, None, None], [0, None, None]]
  star.lam_idx = 0

  # Mask for current chunk 
  m = star.get_masked_chunk(b)

  # The mean training and validation set scatter
  med_training = np.zeros_like(star.lambda_arr)
  med_validation = np.zeros_like(star.lambda_arr)
  
  # Mask transits and outliers
  time = star.time[m]
  flux = star.fraw[m]
  ferr = star.fraw_err[m]
  med = np.nanmedian(flux)
  
  # The precision in the validation set
  validation = [[] for k, _ in enumerate(star.lambda_arr)]

  # The precision in the training set
  training = [[] for k, _ in enumerate(star.lambda_arr)]

  # Setup the GP
  _, amp, tau = star.kernel_params
  gp = george.GP(amp ** 2 * george.kernels.Matern32Kernel(tau ** 2))
  gp.compute(time, ferr)

  # The masks
  masks = list(Chunks(np.arange(0, len(time)), len(time) // star.cdivs))

  # Loop over the different masks
  for i, mask in enumerate(masks):

    print("Section %d/%d..." % (i + 1, len(masks)))

    # Pre-compute (training set)
    pre_t = star.cv_precompute([], b)

    # Pre-compute (validation set)
    pre_v = star.cv_precompute(mask, b)

    # Iterate over lambda
    for k, lam in enumerate(star.lambda_arr):

      # Update the lambda matrix
      star.lam[b][star.lam_idx] = lam

      # Training set
      model = star.cv_compute(b, *pre_t)
      training[k].append(star.fobj(flux - model, med, time, gp, mask))
    
      # Validation set
      model = star.cv_compute(b, *pre_v)
      validation[k].append(star.fobj(flux - model, med, time, gp, mask))

  # Finalize
  training = np.array(training)
  validation = np.array(validation)
  for k, _ in enumerate(star.lambda_arr):

    # Take the mean
    med_validation[k] = np.nanmean(validation[k])
    med_training[k] = np.nanmean(training[k])

  # Compute best model
  i = star.optimize_lambda(validation)
  v_best = med_validation[i]
  t_best = med_training[i]
  star.cdppv_arr[b] = v_best / t_best
  star.lam[b][star.lam_idx] = star.lambda_arr[i]
  print("Found optimum solution at log(lambda) = %.1f." % np.log10(star.lam[b][star.lam_idx]))
  
  lambda_arr = np.array(star.lambda_arr)
  lam = np.array(star.lam)
  lam_idx = star.lam_idx
  
  # Save it!
  np.savez('crossval.npz', lambda_arr = lambda_arr, masks = masks, validation = validation,
           med_training = med_training, med_validation = med_validation, lam = lam,
           lam_idx = lam_idx, v_best = v_best, training = training)

else:
  
  # Load!
  data = np.load('crossval.npz')  
  lambda_arr = data['lambda_arr']
  masks = data['masks']
  validation = data['validation']
  med_training = data['med_training']
  med_validation = data['med_validation']
  lam = data['lam']
  lam_idx = data['lam_idx']
  v_best = data['v_best']
  training = data['training']
  b = 0
     
# Plot!
fig, ax = pl.subplots(4, 1, figsize = (5, 8), sharex = True, sharey = True)
fig.subplots_adjust(hspace = 0.075)

# First x tick will be -infty
lambda_arr[0] = 10 ** (np.log10(lambda_arr[1]) - 3)

# Plot cross-val
for n in range(len(masks)):
  ax[n].plot(np.log10(lambda_arr), validation[:,n], 'r-', alpha = 1)
  ax[n].plot(np.log10(lambda_arr), training[:,n], 'b-', alpha = 1)
  minv = np.argmin(validation[:,n])
  ax[n].annotate('', xy = (np.log10(lambda_arr[minv]), validation[:,n][minv]),
                xycoords = 'data', xytext = (0, 25), textcoords = 'offset points',
                arrowprops = dict(arrowstyle = "-|>", color = 'r'))

# Plot mean
ax[-1].plot(np.log10(lambda_arr), med_training, 'b-', lw = 1., alpha = 1)
ax[-1].plot(np.log10(lambda_arr), med_validation, 'r-', lw = 1., alpha = 1)            
minv = np.argmin(med_validation)
ax[-1].annotate('', xy = (np.log10(lambda_arr[minv]), med_validation[minv]),
                xycoords = 'data', xytext = (0, 25), textcoords = 'offset points',
                arrowprops = dict(arrowstyle = "-|>", color = 'r'))

# Plot best
for axis in ax:
  axis.axvline(np.log10(lam[b][lam_idx]), color = 'k', ls = '--')

# Labels
ax[0].set_ylabel(r'Scatter 1', fontsize = 16, labelpad = 15)
ax[1].set_ylabel(r'Scatter 2', fontsize = 16, labelpad = 15)
ax[2].set_ylabel(r'Scatter 3', fontsize = 16, labelpad = 15)
ax[3].set_ylabel(r'Mean Scatter', fontsize = 16, labelpad = 15)

# Axis range
hi = np.max(validation[0])
lo = np.min(training)
rng = (hi - lo)
for axis in ax:
  axis.set_ylim(lo - 0.15 * rng, hi + 0.15 * rng)
  
# Fix the x ticks
xticks = [np.log10(lambda_arr[0])] + list(np.linspace(np.log10(lambda_arr[1]), np.log10(lambda_arr[-1]), 8))
for axis in ax:
  axis.set_xticks(xticks)
  axis.set_xticklabels(['' for x in xticks])
  axis.set_xlim(np.log10(lambda_arr[0]), 12)

# A hack to mark the first xtick as -infty
labels = ['%.1f' % x for x in xticks]
labels[0] = r'$-\infty$'
ax[-1].set_xticklabels(labels) 
ax[-1].set_xlabel(r'Log $\lambda_1$', fontsize = 16, labelpad = 15)

# Sizes
for axis in ax:
  for tick in axis.get_xticklabels() + axis.get_yticklabels():
    tick.set_fontsize(16)

# Legend
ax[0].plot([-1,-2], [-1,-2], 'b-', label = 'Training')
ax[0].plot([-1,-2], [-1,-2], 'r-', label = 'Validation')
ax[0].legend(loc = 'upper left', fontsize = 10)

# Save
fig.savefig('crossval.pdf', bbox_inches = 0)