#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
injections.py
-------------

'''

from everest.missions.k2 import InjectionStatistics

fig, ax = InjectionStatistics(6.0, show = False)
fig.savefig('injections.pdf', bbox_inches = 'tight')