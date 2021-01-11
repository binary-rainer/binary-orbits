#!/usr/bin/env python
# coding=latin1
"""
TripleFit.py
Created:     Sun Jan 10 18:04:00 2021 by Koehler@Rita
Last change: Mon Jan 11 13:24:59 2021

Python-script to fit Tribble-Orbit
"""

import time
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from astropy.time import Time
from astropy.io  import fits
from KeplerOrbit import KeplerOrbit
from read_tribble import read_tribble


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

    def __init__(self, obs):
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

        return self.compute_chisq()


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

        iPer = 0
        for period in np.logspace( np.log10(P_start), np.log10(P_end), P_npnt, endpoint=False):
            print("========================================================================")
            #print("=== Period",iPer,period)

            if iPer > 0:
                print(iPer*100/P_npnt,"% done, finish in ",end='')
                t= (time.perf_counter() - time0) * (P_npnt-iPer)/iPer/60
                if t < 90: print("{0:.2f} minutes".format(t))
                else:      print("{0:.1f} hours".format(t/60))

            iEcc = 0
            for eccent in np.linspace(e_start, e_end, e_npnt, endpoint=False):
                #print("   === Eccentricity",iEcc,eccent)

                iRat = 0
                for ratio in np.linspace(f_start, f_end, f_npnt, endpoint=False):

                    # create a new orbit-object and put it into The Matrix
                    self.orbit = KeplerOrbit(T0_start,period,1.,eccent,0.,0.,0.)
                    self.orbit.ratio = ratio	# add fr.mass ratio to KeplerOrbit
                    self.orbit.chisq = 1e38	# add chisq-field to KeplerOrbit, set very high
                    self.matrix[iPer,iEcc,iRat] = self.orbit

                    self.xyCM = self.xyOut + ratio*self.xyIn
                    self.exyCM= np.sqrt( self.exyOut*self.exyOut + ratio*self.exyIn*ratio*self.exyIn)


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

                            if chisq < self.matrix[iPer,iEcc,iRat].chisq:
                                self.orbit.chisq = chisq
                                self.matrix[iPer,iEcc,iRat] = self.orbit
                                self.chisq_matrix[iPer,iEcc,iRat]= chisq
                                T0_cent= T0

                                if chisq < chisqmin:
                                    # better than global minimum!
                                    #      T0      P       a       e       o       O       i
                                    print("Better: ",end="")
                                    self.orbit.print_me()
                                    print("	f {0:4.2f}, chi^2 {1:10.1f}".format(ratio,chisq))
                                    chisqmin= chisq
                                    best_iPer,best_iEcc,best_iRat = iPer,iEcc,iRat
                                    print("Index of best orbit:", best_iPer, best_iEcc, best_iRat)

                                # endif chisq<chisqmin
			# endfor T0
                    # endwhile T0_step
                    #print("Final T0-step:",T0_step)

                    iRat += 1
                    # endfor ratio

                if iEcc % 10 == 0:
                    self.orbit.print_me()
                    print(" {0:>10.1f}".format(self.orbit.chisq))

                iEcc += 1
	    # endfor eccent

            iPer += 1
        # endfor period

        print("svdfit finished")
        self.orbit = self.matrix[best_iPer,best_iEcc,best_iRat]
        print("Index of best orbit:", best_iPer, best_iEcc, best_iRat)
        print("Best orbit found: ",end="")
        self.orbit.print_me()
        print("	f {0:4.2f}, chi^2 {1:10.1f}".format(self.orbit.ratio,chisq))
        ratios = np.linspace(f_start, f_end, f_npnt, endpoint=False)
        print("Best ratio from index:", ratios[best_iRat])

        t= time.perf_counter() - time0
        print("Time:",t/60.," minutes")
        plt.ioff()
        return self.orbit

    # end svdfit

    ##################################################################

    def write_svdfit_result(self,filename):
        """write result of svdfit into fits file"""

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

        hdu = fits.PrimaryHDU(data)
        hdu.writeto(filename,overwrite=True)

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
