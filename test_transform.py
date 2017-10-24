#!/usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np
from numpy import sin, cos, arctan2, arcsin
import matplotlib.pyplot as p
from astropy import units as u

l0, b0 = 280.47*u.deg, -32.75*u.deg

def gal2mag(l, b, lp=188.5*u.deg, bp=-7.5*u.deg, anode=False):

    l = l*u.deg
    b = b*u.deg

    Y = -cos(b0)*sin(l0-lp)
    X = sin(b0)*cos(bp) - cos(b0)*sin(bp)*cos(l0-lp)
    Lambda = 90*u.deg + np.arctan2(Y,X)
    Beta = np.arcsin(sin(b0)*sin(bp) + cos(b0)*cos(bp)*cos(l0-lp))
    An = 360*u.deg - Lambda

    Y = -cos(b)*sin(l-lp)
    X = sin(b)*cos(bp) - cos(b)*sin(bp)*cos(l-lp)
    Lambda = 90*u.deg + np.arctan2(Y,X)
    Beta = np.arcsin(sin(b)*sin(bp) + cos(b)*cos(bp)*cos(l-lp))

    Lambda = Lambda + An
    if anode:
        return An
    else:
        return Lambda, Beta

def mag2gal(l, b, lp=188.5*u.deg, bp=-7.5*u.deg, anode=gal2mag(0,0,anode=True)):

    l = l*u.deg
    b = b*u.deg
    l = l - anode

    lpp = 90*u.deg
    bpp = 90*u.deg

    Y = -cos(b)*sin(l-lpp)
    X = sin(b)*cos(bp) - cos(b)*sin(bp)*cos(l-lpp)
    Lambda = lp + np.arctan2(Y,X)
    Beta = np.arcsin(sin(b)*sin(bp) + cos(b)*cos(bp)*cos(l-lpp))
    return Lambda, Beta

if __name__=='__main__':

    from astropy import coordinates as coord
    import MagellanicCoordinates
    from MagellanicCoordinates import Magellanic

    L = np.linspace(0,360,360)
    B = np.zeros_like(L)
    l,b = mag2gal(L,B)
    b = np.rad2deg(b)

    from astropy.coordinates import SkyCoord
    coo = SkyCoord(l,b, frame='galactic')
    myra, mydec = coo.icrs.ra.deg, coo.icrs.dec.deg

    mag = Magellanic(Lambda=0*u.deg, Beta=0*u.deg)
    icrs = mag.transform_to(coord.ICRS)
    gal = mag.transform_to(coord.Galactic)

    gall0, galb0 = gal.l.deg, gal.b.deg

    mag = Magellanic(Lambda=np.linspace(0, 2*np.pi, 360)*u.radian,
                      Beta=np.zeros(360)*u.radian)
    gal = mag.transform_to(coord.Galactic)
    icrs = mag.transform_to(coord.ICRS)

    ra, dec = icrs.ra.deg, icrs.dec.deg

    p.plot(ra, dec, 'k.')
    p.plot(myra, mydec, 'r.')

    p.show()






