#!/usr/bin/env python
# coding=latin1
"""
TripleFit.py
Created:     Sun Jan 10 18:04:00 2021 by Koehler@Rita
Last change: Fri Jan 15 21:53:16 2021

Python-script to fit Tribble-Orbit

2021-Jan-15: fixed bug in find_best_T0() that skipped loop for large steps
"""

import time
import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner

from datetime import datetime
from astropy.time import Time
from astropy.io  import fits
from KeplerOrbit import KeplerOrbit
from read_obs import read_obs

##################################################################

def read_tribble(fname_inner, fname_outer, verbose=False):

    (nameI,obsI) = read_obs(fname_inner, verbose=verbose)
    (nameO,obsO) = read_obs(fname_outer, verbose=verbose)

    cnt = 0

    for iIn in range(len(obsI)):
        if verbose:
            print(iIn,obsI[iIn])

        iOut = (np.where( obsI["MJD"][iIn] == obsO["MJD"]))[0]
        if iOut.size <= 0:
            if verbose:
                print("no match, skipping point")
        else:
            iOut = iOut[0]
            if verbose:
                print(">>>", iOut, obsO[iOut])

            In = obsI[iIn]
            Out = obsO[iOut]

            this = np.array((In["Ref"], In["MJD"],
                             In["Sep"], In["eSep"], In["PA"], In["ePA"],
                             Out["Sep"], Out["eSep"], Out["PA"], Out["ePA"] ),
                            dtype=[("Ref", np.str, 32),
                                   ("MJD", 'f8'),
                                   ('SepI', 'f4'), ('eSepI', 'f4'), ('PA_I', 'f4'), ("ePA_I", 'f4'),
                                   ('SepO', 'f4'), ('eSepO', 'f4'), ('PA_O', 'f4'), ("ePA_O", 'f4')])
            if cnt:
                obs = np.append(obs, this)
            else:
                obs = this
            cnt = cnt+1

    return(nameI, nameO, obs)

##################################################################
# passing only one arg to ensembleSampler doesn't work...
# missing: priors for orbital elements

def log_probability( el, orbf, dummy):
    #print("Elements:",el)
    PeriMJD, period, axis, eccent, periast, node, inclin, ratio = el

    # Eccentricity out of range
    if eccent < 0 or eccent >= 1:
        return -np.inf

    orbf.orbit = KeplerOrbit( PeriMJD, period, axis, eccent, periast, node, inclin)
    orbf.xyCM  = orbf.xyOut + ratio*orbf.xyIn
    orbf.exyCM = np.sqrt( orbf.exyOut*orbf.exyOut + ratio*orbf.exyIn*ratio*orbf.exyIn)

    chisq = orbf.compute_chisq()

    return -0.5 * chisq 	# ln(p)


##################################################################
# we should test this...

def xyErr(sep,eSep, PA,ePA):
    """
    compute error of x/y from Sep,PA
    PA,ePA in Radian
    """

    sinPA = np.sin(PA)
    cosPA = np.cos(PA)

    dx1 = eSep * sinPA
    dy1 = eSep * cosPA
    dx2 = ePA * sep * cosPA
    dy2 = ePA * sep * sinPA

    return np.sqrt( np.array( [ dx1*dx1 + dx2*dx2, dy1*dy1 + dy2*dy2 ] ))


##################################################################

class TripleFit:
    """fit Kepler orbits"""

    def __init__(self, obs, Mass, eMass):
        """
        initialize triple orbit fit from observations
        Mass,eMass in units of axis^3/period^2
        """

        self.Ref = obs["Ref"]
        self.MJD = obs["MJD"]
        self.SepI = obs["SepI"]
        self.eSepI= obs["eSepI"]
        self.PA_I = obs["PA_I"]
        self.ePA_I= obs["ePA_I"]
        self.SepO = obs["SepO"]
        self.eSepO= obs["eSepO"]
        self.PA_O = obs["PA_O"]
        self.ePA_O= obs["ePA_O"]

        self.PA_I[ self.PA_I<0 ] += 360.	# shift into range [0,360]
        self.PA_O[ self.PA_O<0 ] += 360.

        self.PAradI = self.PA_I * np.pi/180.
        self.PAradO = self.PA_O * np.pi/180.

        self.xyIn = np.array([ self.SepI * np.sin(self.PAradI), self.SepI * np.cos(self.PAradI) ])
        self.xyOut= np.array([ self.SepO * np.sin(self.PAradO), self.SepO * np.cos(self.PAradO) ])

        self.exyIn = xyErr( self.SepI, self.eSepI, self.PAradI, self.ePA_I * np.pi/180.)
        self.exyOut= xyErr( self.SepO, self.eSepO, self.PAradO, self.ePA_O * np.pi/180.)

        self.Mass = Mass
        self.eMass = eMass

        self.starname= ""
        self.orbit= None
        self.ratio= None


    ##################################################################

    def compute_chisq(self):
        """compute chi-squared of current orbit"""

        xexp,yexp = self.orbit.MJD_to_xy( self.MJD)
        xyexp = np.array( [xexp,yexp] )		# can we combine these 2 lines?

        dxy = (self.xyCM-xyexp) / self.exyCM

        chisq = np.sum( dxy*dxy )

        if self.Mass:
            model_M = self.orbit.Axis*self.orbit.Axis*self.orbit.Axis / (self.orbit.Period*self.orbit.Period)
            dm = (model_M - self.Mass) / self.eMass
            chisq += dm*dm

        return chisq


    ##################################################################

    def solve_svd(self):

        EAnom = self.orbit.compute_EAnom( self.MJD)

        Amatrix= np.array( [ np.cos(EAnom) - self.orbit.Eccent,
                             np.sqrt(1.-self.orbit.Eccent*self.orbit.Eccent) * np.sin(EAnom) ])
	# Amatrix has shape Nobs,2

        errMat = np.array( [self.exyCM[0,:], self.exyCM[0,:]] )	# same shape as Amatrix
        Xti = np.linalg.lstsq( np.transpose( Amatrix/errMat), self.xyCM[0,:]/self.exyCM[0,:])[0]
        self.orbit.X0, self.orbit.X1 = Xti

        errMat = np.array( [self.exyCM[1,:], self.exyCM[1,:]] )
        Yti = np.linalg.lstsq( np.transpose( Amatrix/errMat), self.xyCM[1,:]/self.exyCM[1,:])[0]
        self.orbit.Y0, self.orbit.Y1 = Yti

        self.orbit.ThieleInnes_to_elements()
        self.orbit.chisq = self.compute_chisq()

        return self.orbit.chisq


    ##################################################################

    def find_best_T0(self, T0_start, T0_npnt, T0_mstp):
        """ find best T0 for position of CM in self.xyCM and elements in self.orbit """

        best_chisq = 1e38	# start very bad

        T0_cent = T0_start
        T0_step = self.orbit.Period*365.25*10./T0_npnt	# will be set to correct step within loop

        while T0_step > T0_mstp:
            T0_step /= 10.
            #print("      T0 step:",T0_step)
            T0_hrng = T0_step*T0_npnt/2
            if T0_step < T0_mstp*2:
                T0_step= T0_mstp	# adjust a little to avoid unnecessary small steps

            for T0 in np.arange( T0_cent-T0_hrng, T0_cent+T0_hrng, T0_step):
                self.orbit.PeriMJD = T0

                chisq = self.solve_svd()

                if chisq < best_chisq:
                    best_chisq = chisq
                    best_orbit = self.orbit
                    #print("new best T0: ",end="")
                    #best_orbit.print_me()
                    #print(" {0:.3f}".format(best_chisq))

        #print("Final T0-step:",T0_step)
        self.orbit = best_orbit


    ##################################################################

    def svdfit(self,
               T0_start, T0_npnt, T0_mstp,
               P_start, P_end, P_npnt,
               e_start, e_end, e_npnt,
               f_start, f_end, f_npnt,
               plot=True):
        """ grid search in T0,p,e,f; SVD for Thiele-Innes """

        print(P_npnt," periods,",e_npnt," eccentricities",f_npnt," fract.mass ratios")
        plt.ion()

        self.matrix= np.empty( (P_npnt,e_npnt,f_npnt), dtype=object)
        self.chisq_matrix = np.zeros( (P_npnt,e_npnt,f_npnt), dtype=np.float64)

        chisqmin = 1e38		# stores global minimum, initially very high
        best_iPer= -1
        best_iEcc= -1
        best_iRat= -1

        time0 = time.perf_counter()

        iRat = 0
        for ratio in np.linspace(f_start, f_end, f_npnt, endpoint=False):
            print("========================================================================")
            if iRat > 0:
                print(iRat*100/f_npnt,"% done, finish in ",end='')
                t= (time.perf_counter() - time0) * (f_npnt-iRat)/iRat/60
                if t < 90: print("{0:.2f} minutes".format(t))
                else:      print("{0:.1f} hours".format(t/60))

            # position of CM according to fract.mass ratio
            self.xyCM = self.xyOut + ratio*self.xyIn
            self.exyCM= np.sqrt( self.exyOut*self.exyOut + ratio*self.exyIn*ratio*self.exyIn)

            iPer = 0
            for period in np.logspace( np.log10(P_start), np.log10(P_end), P_npnt, endpoint=False):
                #print("=== Period",iPer,period)

                iEcc = 0
                for eccent in np.linspace(e_start, e_end, e_npnt, endpoint=False):
                    #print("   === Eccentricity",iEcc,eccent)

                    # create a new orbit-object and put it into The Matrix
                    self.orbit = KeplerOrbit(T0_start,period,1.,eccent,0.,0.,0.)
                    self.orbit.ratio = ratio	# add fract.mass ratio to KeplerOrbit

                    self.find_best_T0( T0_start, T0_npnt, T0_mstp)	# this sets self.orbit to the orbit with the best T0

                    self.matrix[iPer,iEcc,iRat] = self.orbit
                    self.chisq_matrix[iPer,iEcc,iRat]= self.orbit.chisq

                    chisq = self.orbit.chisq
                    if chisq < chisqmin:
                        # better than global minimum!
                        #      T0      P       a       e       o       O       i
                        print("Better: ",end="")
                        self.orbit.print_me()
                        print("	f {0:4.2f}, chi^2 {1:10.1f}".format(ratio,chisq))
                        chisqmin= chisq
                        best_iPer,best_iEcc,best_iRat = iPer,iEcc,iRat
                        #print("Index of best orbit:", best_iPer, best_iEcc, best_iRat)

                        # endif chisq<chisqmin

                    if iEcc % 10 == 0:
                        self.orbit.print_me()
                        print("	f {0:4.2f}, chi^2 {1:10.1f}".format(ratio,chisq))

                    iEcc += 1
                    # endfor eccent

                iPer += 1
                # endfor period

            iRat += 1
            # endfor ratio

        print("svdfit finished")
        self.orbit = self.matrix[best_iPer,best_iEcc,best_iRat]
        print("Index of best orbit:", best_iPer, best_iEcc, best_iRat)
        print("Best orbit found: ",end="")
        self.orbit.print_me()
        print("	f {0:4.2f}, chi^2 {1:10.1f}".format(self.orbit.ratio, self.orbit.chisq))

        ratios = np.linspace(f_start, f_end, f_npnt, endpoint=False)
        print("Best ratio from index:", ratios[best_iRat])

        t= time.perf_counter() - time0
        print("Time:",t/60.," minutes")
        plt.ioff()
        return self.orbit

    # end svdfit

    ##################################################################

    def svdfit_results_to_array(self):
        """ turn array of objects into array of numbers """

        # TODO: test if matrix exists

        shape= self.matrix.shape
        data = np.empty( (9,shape[2],shape[1],shape[0]))
        # 7 orbital elements + ratio + chisq
        # NOTE: pyfits thinks the order of axes should be NAXIS3, -2, -1

        for iPer in range(shape[0]):
            for iEcc in range(shape[1]):
                for iRat in range(shape[2]):
                    el = self.matrix[iPer,iEcc,iRat]

                    #print( type(self.matrix))
                    #print( type(el))
                    #print( el.__class__ )

                    data[0,iRat,iEcc,iPer] = el.PeriMJD
                    data[1,iRat,iEcc,iPer] = el.Period
                    data[2,iRat,iEcc,iPer] = el.Axis
                    data[3,iRat,iEcc,iPer] = el.Eccent
                    data[4,iRat,iEcc,iPer] = el.Periast
                    data[5,iRat,iEcc,iPer] = el.Node
                    data[6,iRat,iEcc,iPer] = el.Inclin
                    data[7,iRat,iEcc,iPer] = el.ratio
                    data[8,iRat,iEcc,iPer] = el.chisq

        return(data)

    #######################################################

    def write_svdfit_result(self,filename):
        """write result of svdfit into fits file"""

        data = self.svdfit_results_to_array()

        hdu = fits.PrimaryHDU(data)
        hdu.writeto(filename,overwrite=True)


    #######################################################

    def read_svdfit_result(self,filename):

        data = fits.getdata(filename)
        shape= data.shape
        print("Data shape:",shape)
	# idx0 are elements, idx1 ratio, idx2 eccent, idx3 period

        self.matrix= np.empty( (shape[3], shape[2], shape[1]), dtype=object)
        self.chisq_matrix = np.zeros( (shape[3], shape[2], shape[1]), dtype=np.float64)

        for iPer in range(shape[3]):
            for iEcc in range(shape[2]):
                for iRat in range(shape[1]):
                    orb = KeplerOrbit( data[0,iRat,iEcc,iPer], #T0
                                       data[1,iRat,iEcc,iPer], #Period
                                       data[2,iRat,iEcc,iPer], #Axis
                                       data[3,iRat,iEcc,iPer], #Eccent
                                       data[4,iRat,iEcc,iPer], #Periast
                                       data[5,iRat,iEcc,iPer], #Node
                                       data[6,iRat,iEcc,iPer] ) #Inclin
                    orb.ratio = data[7,iRat,iEcc,iPer]
                    orb.chisq = data[8,iRat,iEcc,iPer]

                    self.matrix[iPer,iEcc,iRat] = orb
                    self.chisq_matrix[iPer,iEcc,iRat] = orb.chisq

    # end of read_svdfit_result

    ##################################################################

    def MCfit(self, nIter=1000, save_chain=None, save_corner=None):

        nWalker = 48
        nDim = 8	# number of elements + mass ratio

        # obs data is in the object

        svdres = self.svdfit_results_to_array()
        print("shape of svd results:", svdres.shape)

        fchisq = self.chisq_matrix.flatten()

        bestidx = (np.argsort(fchisq))[0:nWalker]	# nWalker indices for best chi-squares

        startpnts = np.zeros( (nWalker, nDim) )

        elnames = [ "T0", "Period", "Axis", "Eccent", "Periast", "Node", "Inclin", "Ratio" ]

        for iEl in range(nDim):
            print("Orbital element #", iEl, elnames[iEl])

            fparam = svdres[iEl,...].flatten()
            startpnts[:,iEl] = fparam[bestidx]

            # shake up T0 and period a little
            if iEl == 0 or iEl == 1:
                startpnts[1:,iEl] += 10. * np.random.randn(nWalker-1)

            print(startpnts[:,iEl])

        print(startpnts.shape)

        sampler = emcee.EnsembleSampler( nWalker, nDim, log_probability, args=(self,0))
        sampler.run_mcmc(startpnts, nIter, progress=True);

        orbel_names = [ "T0", "period", "axis", "e", "omega", "Omega", "i", "ratio" ]
        samps = sampler.get_chain(flat=True)
        #print(samps.shape)

        if save_chain:
            hdu = fits.PrimaryHDU(samps)
            hdu.writeto( save_chain, overwrite=True)

        for iEl in range(8):
            quant = np.percentile(samps[:,iEl], [16,50,84])
            di = np.diff(quant)
            print("{0:>6}: {1:8.2f} +{2:8.2f} -{3:8.2f}".format(orbel_names[iEl], quant[1], di[1], di[0]))

        print("Period:")
        plt.plot(samps[:,1], ".")	# 1=period
        plt.ylabel("Period [years]")
        plt.show()

        fig = plt.figure( figsize=(12,12))
        fig = corner.corner(samps, fig=fig, labels=orbel_names)
        if save_corner:
            plt.savefig(save_corner)

        plt.show()

        return(sampler)


##################################################################

if __name__ == "__main__":

    (nameI,nameO,obs) = read_tribble("TTau-SaSb.pos", "TTau-NSa.pos")

    print("Inner:",nameI)
    print("Outer:",nameO)
    print("Tribble Table:")
    print(obs)

    tribble = TripleFit(obs)
    tribble.svdfit( 50000,100,3., 1e3,1e5,10,  0.0,1.0,10, 0.0, 0.5, 10 )

    tribble.write_svdfit_result("Test.fits")
