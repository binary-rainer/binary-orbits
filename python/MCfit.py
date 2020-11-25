#!/usr/bin/env python
# coding=latin1
"""
MCfit.py
Created:     Tue Nov  3 19:59:55 2020 by Koehler@Rita
Last change: Wed Nov  4 12:45:38 2020

Python-script to do something
"""

import sys
import re
import os.path
import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner

from datetime import datetime
from astropy.time import Time
from read_points import read_points
from KeplerOrbit import KeplerOrbit

#############################################################################
# passing only obs as args to ensembleSampler doesn't work...

def log_probability( el, obs, dummy):
    #print("Elements:",el)
    PeriMJD, period, axis, eccent, periast, node, inclin = el

    # Eccentricity out of range
    if eccent < 0 or eccent >= 1:
        return -np.inf

    orbit = KeplerOrbit( PeriMJD, period, axis, eccent, periast, node, inclin)

    xexp,yexp = orbit.MJD_to_xy(obs["MJD"])

    sepexp = np.sqrt(xexp*xexp + yexp*yexp)
    PA_exp = np.arctan2(xexp,yexp) * 180./np.pi
    PA_exp[ PA_exp<0 ] += 360.	# shift into range [0,360]

    dr = (obs["Sep"]-sepexp) / obs["eSep"]
    dp = (obs["PA"] -PA_exp) / obs["ePA"]

    chisq = np.sum( dr*dr + dp*dp )

    return -0.5 * chisq 	# ln(p) - why?


#############################################################################

for f in sys.argv[1:]:
    print(f)
    starname, obs = read_points(f)

    print(starname)
    print(obs)

    nWalker = 42
    nDim = 7

    startpos = np.array([55348., 6.031, 68.5, 0.79, 255.3, 36.9, 97.0]) + 0.1 * np.random.randn( nWalker, nDim)

    sampler = emcee.EnsembleSampler( nWalker, nDim, log_probability, args=(obs,0))
    sampler.run_mcmc(startpos, 1000, progress=True);

    samps = sampler.get_chain(flat=True)
    print(samps.shape)

    plt.plot(samps[:,2], ".")	# axis
    plt.show()

    labels = [ "T0", "period", "axis", "e", "omega", "Omega", "i" ]
    fig = plt.figure()
    fig = corner.corner(samps, fig=fig, labels=labels)
    plt.show()

#############################################################################
