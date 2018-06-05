# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------------------

import numpy as np
import warnings

from astropy.coordinates.sky_coordinate import SkyCoord
from astropy import units as u
from astropy.time import Time

import matplotlib as mpl

from datetime import datetime, timedelta

'''
 fcmaker: a Python module to automatically create finding charts for ESO OBs in p2.\n
 Copyright (C) 2017,  F.P.A. Vogt
 --- oOo ---
 This file contains general tools for the fcmaker routines.
 Created October 2017, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''

# ----------------------------------------------------------------------------------------
# Design a crosshair custom marker that is sadly missing from maplotlib
def crosshair(inner_r=1, pa = 0):
   ''' 
   The PATH of acrosshair, useful for indicating targets without crowding the field.
   
   Args:
      inner_r: the (relative) inner radius. Default=1
      pa: the position anglem in degrees
   Returns:
      A mpl.path.Path() object.
   '''
   
   verts = [(-1.5, 0),
            (-0.5*inner_r, 0),
            (0, 0.5*inner_r),
            (0, 1.5),
            (0.5*inner_r, 0),
            (1.5, 0),
            (0, -0.5*inner_r),
            (0, -1.5),
            (-1.5, 0),
            (-1.5, 0),
           ]
   
   pa = np.radians(pa)
   rot_mat = np.matrix([[np.cos(pa),-np.sin(pa)],[np.sin(pa),np.cos(pa)]])
   
   for (v, vert) in enumerate(verts):
      verts[v] = (vert*rot_mat).A[0] 
      
   codes = [mpl.path.Path.MOVETO,
            mpl.path.Path.LINETO,
            mpl.path.Path.MOVETO,
            mpl.path.Path.LINETO,
            mpl.path.Path.MOVETO,
            mpl.path.Path.LINETO,
            mpl.path.Path.MOVETO,
            mpl.path.Path.LINETO,
            mpl.path.Path.MOVETO,
            mpl.path.Path.CLOSEPOLY
            ]

   path = mpl.path.Path(verts, codes)
    
   return path

# ----------------------------------------------------------------------------------------
def blind_offset(coord_start, coord_end):
   """
   Compute the offsets between two astropy.sky_coordinate.SkyCoord positions. 
   This function assumes that the sky is flat, and only really makes sense for positions 
   within 15 arcmin.
   
   Args:
      coord_start: a SkyCoord object
      coord_end: a SkyCoord object
      
   Returns:
      The offsets in arcsec
   """

   delta_dec = coord_end.dec-coord_start.dec
   delta_ra = (coord_end.ra-coord_start.ra)*np.cos(start.dec.radian)
   
   if (delta_ra.arcsec > 900*u.arcsec) or (delta_dec.arcsec > 900*u.arcsec):
      warnings.warn('blind_offset: with offsets that large, the sky is no flat anymore!')
   
   return (delta_ra.arcsec, delta_dec.arcsec)

# ----------------------------------------------------------------------------------------
def offset_coord(coord_start, delta_ra=0*u.arcsec, delta_dec=0*u.arcsec):
   """
   Compute a new SkyCoord entity given a starting location and offset. This function 
   assumes that the sky is flat!
   
   Args:
      coord_start: a SKyCoord object
      delta_ra: an RA offset, in u.arcsec
      delta_dec: an Dec offset, in u.arcsec
   
   Returns: 
      a SkyCoord object
   """
   
   if (delta_ra > 900*u.arcsec) or (delta_dec > 900*u.arcsec):
      warnings.warn('offset_coord: with offsets that large, the sky is no flat anymore!')
   
   return SkyCoord(ra = coord_start.ra.deg + delta_ra.to(u.degree).value/np.cos(coord_start.dec.radian),
                   dec = coord_start.dec.deg + delta_dec.to(u.degree).value,
                   unit = (u.degree,u.degree), 
                   frame = 'icrs', 
                   obstime = coord_start.obstime,
                   equinox = coord_start.equinox,
                   )

# ----------------------------------------------------------------------------------------
# ------------ Obsolete after the advent of astropy v3.0 ! -------------------------------
# ----------------------------------------------------------------------------------------
#def propagate_pm(skycoord, epoch, pmRA, pmDec):
#   '''
#   Propagate the proper motion of a target given its epoch and obstime.
#   
#   Args:
#      skycoord (SkyCoord): a SkyCoord item, that must include "obstime"
#      epoch (float): the epoch of the SkyCoord coordinate 
#      pmRa, pmDec (floats): the proper motions in arcsec/year
#   
#   Returns:
#      skycoord: the updated skycoord entry
#      
#   Todo:
#      Implement a proper way of dealing with proper motion, when "velocites" are compatible
#      with SkyCoord in astropy v3.0
#   '''
#   
#   # Convert the epoch (a float) into a datetime entry
#   # https://stackoverflow.com/questions/20911015/decimal-years-to-datetime-in-python
#   epoch_year = int(epoch)
#   rem = epoch - epoch_year
#   
#   base = datetime(epoch_year, 1, 1)
#   epoch = base + rem * timedelta(seconds=(base.replace(year=base.year + 1) - base).total_seconds())
#   
#   # Then compute the time difference bewteen the epoch and the obstime, in year
#   time_dt = skycoord.obstime.to_datetime() - epoch
#   time_dt = time_dt.total_seconds()/u.year.to(u.s)
#   
#   # Now propagate this into the coordinates
#   # This is still quick and dirty, and assumes that the sky is flat ... sigh! 
#   # Maybe better with astropy v3.0 and SkyCoord ???
#   
#   new_skycoord = offset_coord(skycoord, 
#                               delta_ra = pmRA * time_dt*u.arcsec, 
#                               delta_dec = pmDec * time_dt*u.arcsec)
#   
#   
#   return new_skycoord
#
# ---------------------------------------------------------------------------------------- 
def myrotmatrix(angle):
   """
   Defines a clockwise rotation matrix given a certain rotation angle. There is probably
   a smarter way to do this ...
   
   Args:
      angle: a rotation angle in degree
   
   Returns:
      a 2D numpy array
   """
   
   rotmatrix = np.zeros((2,2))
   rotmatrix[0,0] =  np.cos(np.radians(angle))
   rotmatrix[0,1] = -np.sin(np.radians(angle))
   rotmatrix[1,0] =  np.sin(np.radians(angle))
   rotmatrix[1,1] =  np.cos(np.radians(angle))
   return rotmatrix

# ---------------------------------------------------------------------------------------- 
# ------------ NOT USED IN FAVOR OR REAL ASTROPLAN ---------------------------------------
# ---------------------------------------------------------------------------------------- 
#def get_parallactic(target, obsdate,location):
#   '''
#   A function to compute the parallactic angle, given the observing time and telescope.
#   
#   Args:
#      obsdate: a dateime object
#      location: an EarthLocation object
#   Return:
#      The parallactic angle, in degrees
#      
#   Note:
#      Inspired by astroplan's parallactic_angle() function.
#       
#   '''
#
#   # Eqn (14.1) of Meeus' Astronomical Algorithms
#   LST = Time(obsdate).sidereal_time('mean', longitude=location.lon)
#   H = (LST - target.ra).radian
#   q = np.arctan(np.sin(H) / (np.tan(self.location.lat.radian) * 
#       np.cos(coordinate.dec.radian) - np.sin(coordinate.dec.radian)*np.cos(H)))*u.rad
#      
#   return Angle(q)
#
#
# ----------------------------------------------------------------------------------------