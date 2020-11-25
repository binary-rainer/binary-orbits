"""
Python-class for handling a Kepler-orbit

Last Change: Wed Oct  7 09:55:45 2020
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

#######################################################

class KeplerOrbit:
    """classical Kepler orbit"""

    def __init__(self, PeriMJD, Period, Axis, Eccent, Periast, Node, Inclin):
        self.PeriMJD= PeriMJD
        self.Period = Period
        self.Axis   = Axis
        self.Eccent = Eccent
        self.Periast= Periast
        self.Node   = Node
        self.Inclin = Inclin
        self.elements_to_ThieleInnes()

    ############################################

    def print_me(self):
        #      T0       P        a        e        o        O       i
        print("{0:8.2f} {1:6.2f} {2:6.2f} {3:4.2f} {4:6.1f} {5:6.1f}{6:6.1f}".format(
                self.PeriMJD, self.Period, self.Axis, self.Eccent,
                self.Periast, self.Node, self.Inclin), end="")


    ############################################

    def elements_to_ThieleInnes(self):
        peri = self.Periast * np.pi/180.
        node = self.Node    * np.pi/180.
        inc  = self.Inclin  * np.pi/180.

        cp = np.cos(peri)
        sp = np.sin(peri)
        cn = np.cos(node)
        sn = np.sin(node)
        ci = np.cos(inc)

        self.X0 = self.Axis * ( cp * sn + sp * cn * ci)
        self.Y0 = self.Axis * ( cp * cn - sp * sn * ci)
        self.X1 = self.Axis * (-sp * sn + cp * cn * ci)
        self.Y1 = self.Axis * (-sp * cn - cp * sn * ci)

    ############################################

    def ThieleInnes_to_elements(self):
        # might be faster to compute sin(o±O) and cos(o±O) here and store them

        opO = np.arctan2( self.X0-self.Y1, self.X1+self.Y0 )	# omega + Omega
        omO = np.arctan2(-self.X0-self.Y1,-self.X1+self.Y0 )	# omega - Omega

        self.Periast= (opO + omO) * 180./2./np.pi	# omega = argument of periastron
        self.Node   = (opO - omO) * 180./2./np.pi	# Omega = P.A. of node

        if self.Periast < 0.: self.Periast += 360.
        if self.Node    < 0.: self.Node += 360.

        # from astrometry, we cannot distinguish the two nodes
        # Hilditch, p.52: convention demands we adopt the solution for which Omega < pi
        if self.Node >= 180.:
            self.Node -= 180.
            self.Periast -= 180.

        # avoid singularity for sin(opO) = 0
        #
        if np.fabs(np.sin(opO)*np.sin(omO)) > np.fabs(np.cos(opO)*np.cos(omO)):
            help  = -((self.X0+self.Y1) * np.sin(opO)) / ((self.X0-self.Y1) * np.sin(omO))
            cosinc= (1.-help) / (1.+help)
            self.Axis= (self.X0-self.Y1) / np.sin(opO) / (1.+cosinc)
        else:
            help  = ((self.X1-self.Y0) * np.cos(opO)) / ((self.X1+self.Y0) * np.cos(omO))
            cosinc= (1.+help) / (1.-help)
            self.Axis= (self.X1+self.Y0) / np.cos(opO) / (1.+cosinc)

        self.Inclin= np.arccos(cosinc) * 180./np.pi


    ############################################

    def compute_EAnom(self,MJD):
        # only elliptic orbits for now

        M = (MJD - self.PeriMJD) / (self.Period*365.25)	# mean anomaly in revolutions
        M,_ = np.modf(M)		# need only fractional part
        M  *= 2.*np.pi			# in radian in [0,2pi]

        if self.Eccent == 0.:
            return M

        if self.Eccent >= 1.:
            print("Sorry, only elliptic orbits have been implemented!")
            return 0.

        i = 0
        E = M
        while i<100:
            i += 1
            E0 = E
            E  = M + self.Eccent * np.sin(E)
            # make sure max sees an array, even if we have a scalar
            if max( np.append(abs(E-E0),0.)) < 1e-7:
                break
        return E

    ############################################

    def rad_from_cnu(self, cnu):
        # cnu = cos(true Anomaly)

        if self.Eccent != 1.:
            return abs(1. - self.Eccent * self.Eccent) / (1. + self.Eccent * cnu)
	    # Montenbruck p. 54/63 and 124, a = 1 with Thiele-Innes
        else:
            return 2. / (1. + cnu)

    ############################################
    # nu in degrees!

    def true_anom_to_xy(self, nu):
        nur = nu * np.pi/180.

        cnu = np.cos(nur)
        snu = np.sin(nur)

        r = self.rad_from_cnu(cnu)

        x = r * (self.X0 * cnu + self.X1 * snu)
        y = r * (self.Y0 * cnu + self.Y1 * snu)

        return x,y

    ############################################
    # result is in the length units of the orbital elements
    # divided by period converted to seconds (e.g. mas/s)
    # based on, e.g., Hilditch p.42

    def true_anom_to_rv(self,nu):

        K = 2.*np.pi * self.Axis * np.sin(self.Inclin*np.pi/180.) \
            / (self.Period*365.25*86400.) / np.sqrt(1.-self.Eccent*self.Eccent)

        om= self.Periast * np.pi/180.	# lowercase omega in rad

        rv = K * (self.Eccent * np.cos(om) + np.cos(nu+om))
        return rv * 1.


    ############################################

    def MJD_to_xy(self, MJD):
        E = self.compute_EAnom(MJD)

        if self.Eccent < 1.:
            rcnu = np.cos(E) - self.Eccent
            rsnu = np.sqrt(1.-self.Eccent*self.Eccent) * np.sin(E)
        elif self.Eccent == 1.:
            # E is in fact true anomaly here
            r    = 2./(1.+np.cos(E))
            rcnu = r * np.cos(E)
            rsnu = r * np.sin(E)
        else:
            rcnu = self.Eccent - np.cosh(E)
            rsnu = np.sqrt(self.Eccent*self.Eccent - 1.) * np.sinh(E)

        x = self.X0 * rcnu + self.X1 * rsnu
        y = self.Y0 * rcnu + self.Y1 * rsnu

        return x,y

    ############################################

    def MJD_to_z(self, MJD):

        E = self.compute_EAnom(MJD)

        if self.Eccent < 1.:
            rcnu = np.cos(E) - self.Eccent
            rsnu = np.sqrt(1.-self.Eccent*self.Eccent) * np.sin(E)
        elif self.Eccent == 1.:
            # E is in fact true anomaly here
            r    = 2./(1.+np.cos(E))
            rcnu = r * np.cos(E)
            rsnu = r * np.sin(E)
        else:
            rcnu = self.Eccent - np.cosh(E)
            rsnu = np.sqrt(self.Eccent*self.Eccent - 1.) * np.sinh(E)

        peri = self.Periast * np.pi/180.
        inc  = self.Inclin  * np.pi/180.

        Z0 = self.Axis * np.sin(peri) * np.sin(inc)
        Z1 = self.Axis * np.cos(peri) * np.sin(inc)	# buf fixed 7-Oct-2020

        z = Z0 * rcnu + Z1 * rsnu

        return z

    ############################################

    def MJD_to_rv(self, MJD):
        E = self.compute_EAnom(MJD)

        if self.Eccent < 1.:
            nu = 2. * np.arctan( np.sqrt( (1.+self.Eccent)/(1.-self.Eccent) ) * np.tan(E/2.))
        elif self.Eccent == 1.:
            nu = E	# E is in fact true anomaly here
        else:
            nu = 2. * np.arctan( np.sqrt( (self.Eccent+1.)/(self.Eccent-1.) ) * np.tanh(E/2.))

        if self.Eccent >= 1.:
            print("Warning: I have no idea whether rv for non-elliptic orbits are correct!")

        return self.true_anom_to_rv(nu)


#############################################################################
