#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
cbv.py
------

'''

import everest
import os
import matplotlib.pyplot as pl
import numpy as np

def GetChunk(time, breakpoints, b, mask = []):
  '''
  Returns the indices corresponding to a given light curve chunk.
  
  :param int b: The index of the chunk to return

  '''

  M = np.delete(np.arange(len(time)), mask, axis = 0)
  if b > 0:
    res = M[(M > breakpoints[b - 1]) & (M <= breakpoints[b])]
  else:
    res = M[M <= breakpoints[b]]
  return res

# Set up the plot
fig, ax = pl.subplots(3, 3, figsize = (13, 10))
ax = ax.flatten()

# Loop over all campaigns
for campaign in range(9):
  
  # Get the data
  t, x1, x2 = np.loadtxt(os.path.join('cbv', '%02d.tsv' % campaign), unpack = True)
  with open(os.path.join('cbv', '%02d.tsv' % campaign)) as f:
    breakpoints = [int(b) for b in f.readline().replace("#", '').split(',')]
  
  if campaign == 0:
    bad = np.arange(664, 1912)
    x1[bad] = np.nan
    x2[bad] = np.nan
  
  # Plot the first two CBVs in each chunk of this quarter
  offset = 0
  for b in range(len(breakpoints)):
    i = GetChunk(t, breakpoints, b)
    ax[campaign].plot(t[i] + offset, (x2[i] - np.nanmin(x2[i])) / (np.nanmax(x2[i]) - np.nanmin(x2[i])), color = 'r', alpha = 0.5)
    ax[campaign].plot(t[i] + offset, (x1[i] - np.nanmin(x1[i])) / (np.nanmax(x1[i]) - np.nanmin(x1[i])), color = 'b')
    offset += 5
  
  # Appearance
  ax[campaign].set_title('C%02d' % campaign, fontsize = 22)
  ax[campaign].margins(0.1, 0.1)
  ax[campaign].set_xticks([])
  ax[campaign].set_yticks([])
  
# Save
fig.savefig('cbv.pdf', bbox_inches = 'tight')