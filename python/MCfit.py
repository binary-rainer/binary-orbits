#!/usr/bin/env python
# coding=latin1
"""
MCfit.py
Created:     Tue Nov  3 19:59:55 2020 by Koehler@Rita
Last change: Thu Jan  7 18:39:39 2021

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
from astropy.io  import fits
from astropy.time import Time
from astropy.table import Table
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

    # column names from text file:
    #dr = (obs["Sep"]-sepexp) / obs["eSep"]
    #dp = (obs["PA"] -PA_exp) / obs["ePA"]

    # column names from fits file (what was I thinking???)
    dr = (obs["Separation"]    -sepexp) / obs["Sep Error"]
    dp = (obs["Position Angle"]-PA_exp) / obs["PA Error"]

    chisq = np.sum( dr*dr + dp*dp )

    return -0.5 * chisq 	# ln(p) - why?


#############################################################################

nWalker = 42
nDim = 7	# actually, number of elements


# this has different column names than inside class OrbitFit!
obs = Table.read("TTau-SaSb_obs.fits")
print("obs table:")
print(obs.info)


svdres = fits.getdata("SaSb-fit01_fit.fits")
print("shape of svd results:", svdres.shape)


fchisq = svdres[7,...].flatten()

# nWalker indices for best chi-squares
bestidx = (np.argsort(fchisq))[0:nWalker]

startpnts = np.zeros( (nWalker, nDim) )

for iEl in range(nDim):
    print(iEl)

    fparam = svdres[iEl,...].flatten()
    startpnts[:,iEl] = fparam[bestidx]


print(startpnts.shape)


sampler = emcee.EnsembleSampler( nWalker, nDim, log_probability, args=(obs,0))
sampler.run_mcmc(startpnts, 1000, progress=True);

samps = sampler.get_chain(flat=True)
print(samps.shape)

plt.plot(samps[:,2], ".")	# axis
plt.show()

labels = [ "T0", "period", "axis", "e", "omega", "Omega", "i" ]
fig = plt.figure()
fig = corner.corner(samps, fig=fig, labels=labels)
plt.show()

#############################################################################
