"""
Python class for fitting Kepler orbits to astrometric data

Last Change: Sat Oct 19 00:00:44 2019
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.io  import fits
from KeplerOrbit import KeplerOrbit
import time

#############################################################################

class OrbitFit:
    """fit Kepler orbits"""

    """
    obs data is passed as structured array
    """

    def __init__(self, obs):
        self.Ref = obs["Ref"]
        self.MJD = obs["MJD"]
        self.Sep = obs["Sep"]
        self.eSep= obs["eSep"]
        self.PA  = obs["PA"]
        self.ePA = obs["ePA"]

        self.PA[ self.PA<0 ] += 360.	# shift into range [0,360]

        self.PArad= self.PA * np.pi/180.
        self.xobs = self.Sep * np.sin(self.PArad)
        self.yobs = self.Sep * np.cos(self.PArad)

        dx1 = self.eSep * np.sin(self.PArad)
        dy1 = self.eSep * np.cos(self.PArad)
        dx2 = self.ePA * np.pi/180. * self.yobs
        dy2 = self.ePA * np.pi/180. * self.xobs
        self.exobs= np.sqrt( dx1*dx1 + dx2*dx2 )
        self.eyobs= np.sqrt( dy1*dy1 + dy2*dy2 )

        self.orbit= None

    # end of __init__


    def save_observations(self,filename):

        hdr = fits.Header()
        hdr['EXTNAME'] = 'Astrometry'

        ff= fits.BinTableHDU.from_columns(
            [fits.Column(name='Ref', format='16A',	array=self.Ref),
             fits.Column(name='MJD', format='D', 	array=self.MJD),
             fits.Column(name='Separation', format='E', array=self.Sep),
             fits.Column(name='Sep Error',  format='E', array=self.eSep),
             fits.Column(name='Position Angle',format='E',array=self.PA, unit='deg'),
             fits.Column(name='PA Error',     format='E',array=self.ePA,unit='deg'),
             fits.Column(name='x',      format='E', array=self.xobs ),
             fits.Column(name='x error',format='E', array=self.exobs),
             fits.Column(name='y',      format='E', array=self.yobs ),
             fits.Column(name='y error',format='E', array=self.eyobs) ],
            header=hdr )
            # would be nice to have units for sep,x,y

        ff.writeto(filename,clobber=True)


    def plot_orbpnt(self, iPnt, color='orange'):
        PA = self.PArad[iPnt]
        daz= self.Sep[iPnt] * self.ePA[iPnt] * np.pi/180.	# delta azimuth
        s  = self.Sep[iPnt]
        ds = self.eSep[iPnt]

        w = np.arange(361) * np.pi/180.
        xx= (s + ds*np.cos(w)) * np.sin(PA) - daz*np.sin(w) * np.cos(PA)
        yy= (s + ds*np.cos(w)) * np.cos(PA) + daz*np.sin(w) * np.sin(PA)

        plt.plot(xx,yy,color=color)
        #print("Point",iPnt,":",np.mean(xx),np.mean(yy))

        if self.orbit:
            xexp,yexp = self.orbit.MJD_to_xy( self.MJD[iPnt])
            xobs = s * np.sin(PA)
            yobs = s * np.cos(PA)
            plt.plot( [xexp,xobs], [yexp,yobs], linestyle="-", color=color)



    def plot_orbit(self, pnt_color='red', orb_color='blue', orb_lw=1):
        # for dark bg: pnt_color='orange', orb_color='lime', orb_lw=2
        #plt.clf() - this destroys the ellipse showing TTauS' circumbinary disk
        #print("plot_orbit: lw=",orb_lw)
        if self.orbit:
            x,y = self.orbit.true_anom_to_xy( np.arange(361))
            plt.plot(x,y, color=orb_color, lw=orb_lw)
            plt.plot( [0.,x[0]], [0.,y[0]], linestyle="-.",color=orb_color, lw=orb_lw)	# periastron
            xnode,ynode = self.orbit.true_anom_to_xy( np.array([0.,180.]) - self.orbit.Periast)
            plt.plot(xnode,ynode, linestyle="--", color=orb_color, lw=orb_lw)		# line of nodes
        else:
            x,y = self.xobs,self.yobs

        plt.plot(self.xobs, self.yobs, '+', mew=1, color=pnt_color)
        for i in range(self.Sep.shape[0]):
            self.plot_orbpnt(i, color=pnt_color)

        #plt.axis('equal')	# isotropic axes
        plt.gca().set_aspect('equal', adjustable='box')	# from stackoverflow
        #if self.orbit:
        #    plt.xlim( np.nanmax(x), np.nanmin(x))	# invert x-axis
        #else:
        # x,y are the orbit or the observations - make a symmetric square plot
        xymax = np.nanmax(np.fabs([x,y]))*1.1
        plt.xlim( xymax,-xymax)		# invert x-axis
        plt.ylim(-xymax, xymax)

        plt.plot( [0],[0], 'x', mew=2, color=pnt_color)
        plt.draw()
        #plt.show(block=block) # not necessary if ion() is used


    def plot_chisq_matrix(self):
        fig = plt.figure(2)
        plt.clf()
        plt.imshow(np.transpose(self.chisq_matrix),interpolation='nearest',origin='lower',cmap='gray')
        plt.xlabel("Period")
        plt.ylabel("Eccentricity")
        plt.draw()


    def compute_chisq(self):
        """compute chi-squared of current orbit"""

        xexp,yexp = self.orbit.MJD_to_xy( self.MJD)

        sepexp = np.sqrt(xexp*xexp + yexp*yexp)
        PA_exp = np.arctan2(xexp,yexp) * 180./np.pi
        PA_exp[ PA_exp<0 ] += 360.	# shift into range [0,360]

        dr = (self.Sep-sepexp) / self.eSep
        dp = (self.PA -PA_exp) / self.ePA

        chisq = np.sum( dr*dr + dp*dp )
        return chisq

    # end of compute_chisq


    def print_residuals(self):
        """ print observed & expected positions and differences """

        xexp,yexp = self.orbit.MJD_to_xy( self.MJD)

        sepexp = np.sqrt(xexp*xexp + yexp*yexp)
        PA_exp = np.arctan2(xexp,yexp) * 180./np.pi
        PA_exp[ PA_exp<0 ] += 360.	# shift into range [0,360]

        dr = (self.Sep-sepexp) / self.eSep
        dp = (self.PA -PA_exp) / self.ePA

        Epoch = 2000.0 + (self.MJD-51544.5) / 365.25

        print("Julian date  Observed____________  Expected____________  Sigma____ Ref")
        for i in range(Epoch.size):
            print("{0:10.3f}:  {1:5.1f} mas {2:6.1f} deg, {3:5.1f} mas {4:6.1f} deg, {5:4.1f} {6:4.1f} {7}".
                  format( Epoch[i], self.Sep[i], self.PA[i], sepexp[i], PA_exp[i], dr[i], dp[i], self.Ref[i]))

        chisq = np.sum( dr*dr + dp*dp )
        print("chi^2:",chisq)
        print("red.chi^2:", chisq / (2*Epoch.size-7))


    # NOTE: loops stop _before_ value *_end is tried

    def svdfit(self, T0_start,T0_npnt,T0_mstp, P_start,P_end,P_npnt, e_start,e_end,e_npnt,plot=True):
        """ grid search in T0,p,e; SVD for Thiele-Innes """

        print(P_npnt," periods,",e_npnt," eccentricities")
        plt.ion()

        self.matrix= np.empty( (P_npnt,e_npnt), dtype=object)
        self.chisq_matrix = np.zeros( (P_npnt,e_npnt), dtype=np.float64)

        chisqmin = 1e38		# stores global minimum, initially very high
        best_iPer= -1
        best_iEcc= -1

        time0 = time.clock()

        iPer = 0
        for period in np.logspace( np.log10(P_start), np.log10(P_end), P_npnt, endpoint=False):
            print("========================================================================")
            #print("=== Period",iPer,period)

            if iPer > 0:
                print(iPer*100/P_npnt,"% done, finish in ",end='')
                t= (time.clock() - time0) * (P_npnt-iPer)/iPer/60
                if t < 90: print("{0:.2f} minutes".format(t))
                else:      print("{0:.1f} hours".format(t/60))

            iEcc = 0
            for eccent in np.linspace(e_start,e_end,e_npnt, endpoint=False):
                #print("   === Eccentricity",iEcc,eccent)

                # create a new orbit-object and put it into The Matrix
                self.orbit = KeplerOrbit(T0_start,period,1.,eccent,0.,0.,0.)
                self.orbit.chisq= 1e38	# add chisq-field to KeplerOrbit, set very high
                self.matrix[iPer,iEcc] = self.orbit

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

                        EAnom = self.orbit.compute_EAnom( self.MJD)

                        Amatrix= np.array( [ np.cos(EAnom) - self.orbit.Eccent,
                                             np.sqrt(1.-self.orbit.Eccent*self.orbit.Eccent) * np.sin(EAnom) ])
                        # Amatrix has shape Nobs,2

                        errMat = np.array( [self.exobs, self.exobs] )	# same shape as Amatrix
                        Xti = np.linalg.lstsq( np.transpose( Amatrix/errMat), self.xobs/self.exobs)[0]
                        self.orbit.X0, self.orbit.X1 = Xti

                        errMat = np.array( [self.eyobs, self.eyobs] )
                        Yti = np.linalg.lstsq( np.transpose( Amatrix/errMat), self.yobs/self.eyobs)[0]
                        self.orbit.Y0, self.orbit.Y1 = Yti

                        self.orbit.ThieleInnes_to_elements()

                        chisq = self.compute_chisq()

                        if chisq < self.matrix[iPer,iEcc].chisq:
                            self.orbit.chisq = chisq
                            self.matrix[iPer,iEcc] = self.orbit
                            self.chisq_matrix[iPer,iEcc]= chisq
                            T0_cent= T0

                            if chisq < chisqmin:
                                # better than global minimum!
                                #      T0      P       a       e       o       O       i
                                print("Better: ",end="")
                                self.orbit.print_me()
                                print(" {0:>10.1f}".format(chisq))
                                chisqmin= chisq
                                best_iPer,best_iEcc = iPer,iEcc
                                if plot:
                                    fig = plt.figure(1)
                                    plt.clf()
                                    self.plot_orbit()
                                    self.plot_chisq_matrix()
                                    plt.pause(0.01)

                        # endif chisq<chisqmin
                    # endfor T0
                # endwhile T0_step
                #print("Final T0-step:",T0_step)

                if iEcc % 10 == 0:
                    self.orbit.print_me()
                    print(" {0:>10.1f}".format(self.orbit.chisq))

                iEcc += 1
            # endfor eccent

            iPer += 1
        # endfor period

        print("svdfit finished")
        self.orbit = self.matrix[best_iPer,best_iEcc]
        print("Best orbit found: ",end="")
        self.orbit.print_me()
        print(" {0:>10.1f}".format(chisq))
        t= time.clock() - time0
        print("Time:",t/60.," minutes")
        plt.ioff()
        return self.orbit

    # end svdfit


    def write_svdfit_result(self,filename):
        """write result of svdfit into fits file"""

        # TODO: test if matrix exists

        shape= self.matrix.shape
        data = np.empty( (8,shape[1],shape[0]))
        # 7 orbital elements + chisq
        # NOTE: pyfits thinks the order of axes should be NAXIS3, -2, -1

        for iPer in range(shape[0]):
            for iEcc in range(shape[1]):
                el = self.matrix[iPer,iEcc]

                #print( type(self.matrix))
                #print( type(el))
                #print( el.__class__ )

                data[0,iEcc,iPer] = el.PeriMJD
                data[1,iEcc,iPer] = el.Period
                data[2,iEcc,iPer] = el.Axis
                data[3,iEcc,iPer] = el.Eccent
                data[4,iEcc,iPer] = el.Periast
                data[5,iEcc,iPer] = el.Node
                data[6,iEcc,iPer] = el.Inclin
                data[7,iEcc,iPer] = el.chisq

        hdu = fits.PrimaryHDU(data)
        hdu.writeto(filename,overwrite=True)


    def read_svdfit_result(self,filename):
        """read result of svdfit from fits file"""

        data = fits.getdata(filename)
        print("Data's shape:",data.shape)
        e_npnt,P_npnt = data.shape[1:3]
        print(P_npnt,e_npnt)

        self.matrix= np.empty( (P_npnt,e_npnt), dtype=object)
        best_chisq= data[7,0,0]+1.	# trigger copy of el at first iteration

        for iPer in range(P_npnt):
            for iEcc in range(e_npnt):
                el = KeplerOrbit( PeriMJD= data[0,iEcc,iPer],
                                  Period = data[1,iEcc,iPer],
                                  Axis   = data[2,iEcc,iPer],
                                  Eccent = data[3,iEcc,iPer],
                                  Periast= data[4,iEcc,iPer],
                                  Node   = data[5,iEcc,iPer],
                                  Inclin = data[6,iEcc,iPer] )
                el.chisq = data[7,iEcc,iPer]
                self.matrix[iPer,iEcc] = el
                if best_chisq > el.chisq:
                    # set self.orbit to best fit
                    self.orbit= el
                    best_chisq= el.chisq

        self.chisq_matrix = data[7,...]


#############################################################################
