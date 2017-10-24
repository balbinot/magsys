#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
==========================================================
Create a new coordinate class (for the Magellanic stream)
==========================================================

-------------------

*Adapted by Eduardo Balbinot (eduardo.balbinot@gmail.com)

*Using the IDL scripts found here*
https://github.com/bsmartforever/wham/tree/master/PRO/Nidever/mag2gal

*Originally by: Adrian Price-Whelan, Erik Tollerud*
http://docs.astropy.org/en/stable/generated/examples/coordinates/plot_sgr-coordinate-frame.html

*License: BSD*

-------------------

"""

##############################################################################
# Make `print` work the same in all versions of Python, set up numpy,
# matplotlib, and use a nicer set of plot parameters:

import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

##############################################################################
# Import the packages necessary for coordinates

from astropy.coordinates import frame_transform_graph
from astropy.coordinates.matrix_utilities import rotation_matrix, matrix_product, matrix_transpose
import astropy.coordinates as coord
import astropy.units as u

##############################################################################
# The first step is to create a new class, which we'll call
# ``Sagittarius`` and make it a subclass of
# `~astropy.coordinates.BaseCoordinateFrame`:

class Magellanic(coord.BaseCoordinateFrame):
    """

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)

    Lambda : `Angle`, optional, must be keyword
        The longitude-like angle corresponding to the Magellanic stream
    Beta : `Angle`, optional, must be keyword
        The latitude-like angle corresponding to the Magellanic stream
    distance : `Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.

    pm_Lambda_cosBeta : :class:`~astropy.units.Quantity`, optional, must be keyword
        The proper motion along the stream in ``Lambda`` (including the
        ``cos(Beta)`` factor) for this object (``pm_Beta`` must also be given).
    pm_Beta : :class:`~astropy.units.Quantity`, optional, must be keyword
        The proper motion in Declination for this object (``pm_ra_cosdec`` must
        also be given).
    radial_velocity : :class:`~astropy.units.Quantity`, optional, must be keyword
        The radial velocity of this object.

    """
    default_representation = coord.SphericalRepresentation
    default_differential = coord.SphericalCosLatDifferential

    frame_specific_representation_info = {
        coord.SphericalRepresentation: [
            coord.RepresentationMapping('lon', 'Lambda'),
            coord.RepresentationMapping('lat', 'Beta'),
            coord.RepresentationMapping('distance', 'distance')],
        coord.SphericalCosLatDifferential: [
            coord.RepresentationMapping('d_lon_coslat', 'pm_Lambda_cosBeta'),
            coord.RepresentationMapping('d_lat', 'pm_Beta'),
            coord.RepresentationMapping('d_distance', 'radial_velocity')],
        coord.SphericalDifferential: [
            coord.RepresentationMapping('d_lon', 'pm_Lambda'),
            coord.RepresentationMapping('d_lat', 'pm_Beta'),
            coord.RepresentationMapping('d_distance', 'radial_velocity')]
    }

    frame_specific_representation_info[coord.UnitSphericalRepresentation] = \
        frame_specific_representation_info[coord.SphericalRepresentation]
    frame_specific_representation_info[coord.UnitSphericalCosLatDifferential] = \
        frame_specific_representation_info[coord.SphericalCosLatDifferential]
    frame_specific_representation_info[coord.UnitSphericalDifferential] = \
        frame_specific_representation_info[coord.SphericalDifferential]

from sys import argv

#Poles at (l,b) = 188.5, -7.5

MAG_PHI = (180-188.5) * u.degree
MAG_THETA = (90+7.5) * u.degree
MAG_PSI = (90+32.8610) * u.degree

# Generate the rotation matrix using the x-convention (see Goldstein)
A = rotation_matrix(MAG_PHI, (0,0,-1))
B = rotation_matrix(MAG_THETA, (0,-1,0))
G = rotation_matrix(MAG_PSI, (0,0,-1))
R = np.diag([1.,1.,-1.])
MAG_MATRIX = matrix_product(R,G,B,A)

@frame_transform_graph.transform(coord.StaticMatrixTransform, coord.Galactic, Magellanic)
def galactic_to_magellanic():
    """ Compute the transformation matrix from Galactic spherical to
        heliocentric Magellanic coordinates.
    """
    return MAG_MATRIX

@frame_transform_graph.transform(coord.StaticMatrixTransform, Magellanic, coord.Galactic)
def magellanic_to_galactic():
    """ Compute the transformation matrix from heliocentric Magellanic
        coordinates to spherical Galactic.
    """
    return matrix_transpose(MAG_MATRIX)


if __name__=='__main__':
    mag = Magellanic(Lambda=0*u.deg, Beta=0*u.deg)
    icrs = mag.transform_to(coord.ICRS)
    gal = mag.transform_to(coord.Galactic)

    gall0, galb0 = gal.l.deg, gal.b.deg

    print "Magellanic (%.4f, %.4f) = ICRS (%.4f, %.4f) = Galactic (%.4f, %.4f)" % (mag.Lambda.deg, mag.Beta.deg, icrs.ra.deg, icrs.dec.deg, gal.l.deg, gal.b.deg)

    mag = Magellanic(Lambda=np.linspace(-2*np.pi, 2*np.pi, 256)*u.radian,
                      Beta=np.zeros(256)*u.radian)
    gal = mag.transform_to(coord.Galactic)
    icrs = mag.transform_to(coord.ICRS)

    fig, axes = plt.subplots(2, 1, figsize=(8, 10))

    axes[0].set_title("Magellanic")
    axes[0].plot(mag.Lambda.wrap_at(180*u.deg).deg, mag.Beta.deg,
                 linestyle='none', marker='.')

    axes[1].set_title("Galactic")
    axes[1].plot(gal.l.deg, gal.b.deg, ls='none',
                 marker='.')

    axes[1].plot(gall0, galb0, 'ro', ms=8)
    axes[1].set_xlim(360,0)
    axes[1].axhline(-32.75)

    plt.show()

