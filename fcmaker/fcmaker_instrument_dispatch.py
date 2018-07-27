# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------------------

import os
import numpy as np
import warnings

from astropy.coordinates.sky_coordinate import SkyCoord
from astropy import units as u
from astropy.time import Time
import copy

from datetime import datetime
from datetime import timedelta
from dateutil import parser as dup

from . import fcmaker_metadata as fcm_m
from . import fcmaker_tools as fcm_t
from . import fcmaker_muse as fcm_muse
from . import fcmaker_hawki as fcm_hawki
from . import fcmaker_xshooter as fcm_xshooter


'''
 fcmaker: a Python module to automatically create finding charts for ESO OBs in p2.\n
 Copyright (C) 2017,  F.P.A. Vogt
 --- oOo ---
 This file contains general functions dispatching to specific instrment functions for the 
 fcmaker routines. 
 Created October 2017, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''

# ----------------------------------------------------------------------------------------
def get_bk_image(fc_params):
   '''
   Returns the default SkyView survey name for all supported instruments.
   
   Args:
      fc_params: the parameters of the finding chart (a dictionnary)
   Returns:
      The survey name as a string
   '''
   
   # Define which default background image parameters to use for specific instruments
   # Use Skyview to download it from well known sources.
   # See http://astroquery.readthedocs.io/en/latest/skyview/skyview.html
   # To see all possible catalogues, type: SkyView.survey_dict['overlay_blue']
   
   if fc_params['inst'] == 'MUSE':
      if 'NFM' in fc_params['ins_mode']:
         return fcm_muse.bk_image_nfm
      else:
         return fcm_muse.bk_image_wfm
      
   elif fc_params['inst'] == 'HAWKI':
      # For HAWKI, in case of J, H or K, then use the corresponding 2MASS image
      if fc_params['acq']['filter'] in ['J','H','Ks']:
         return '2MASS-'+fc_params['acq']['filter'][0]
      else:
         return fcm_hawki.bk_image
         
   elif fc_params['inst'] =='XSHOOTER':
      return fcm_xshooter.bk_image
   
   else:
      raise Exception('Ouch! Instrument unknown ...') 

# ----------------------------------------------------------------------------------------
def get_fields_dict(fc_params):
   '''
   Returns the dictionnary of the field plotting parameters, given certain fc parameters.
   
   Args:
      fc_params: the parameters of the finding chart (a dictionnary)
   Returns:
      A dictionary of lists.
   '''
   
   if fc_params['inst'] == 'MUSE':
      return fcm_muse.get_fields_dict(fc_params)
   
   elif fc_params['inst'] == 'HAWKI':
      return fcm_hawki.get_fields_dict(fc_params)
   
   elif fc_params['inst'] == 'XSHOOTER':
      return fcm_xshooter.get_fields_dict(fc_params)
      
   else:
      raise Exception('Ouch! Instrument unknown ...')   

# ----------------------------------------------------------------------------------------
def plot_field(ax1, ax2, fc_params, field):
   '''
   Triggers the plotting of a specific field.
   
   Args:
      ax1,ax2: the left and right plots 'axis'
      fc_params: the dictionnary of parameters for the OB
      field: the specific field list of parameters
   '''
   
   if fc_params['inst'] == 'MUSE':
      fcm_muse.plot_field(ax1,ax2,fc_params,field)
   
   elif fc_params['inst'] == 'HAWKI':
      fcm_hawki.plot_field(ax1,ax2,fc_params,field)
      
   elif fc_params['inst'] == 'XSHOOTER':
      fcm_xshooter.plot_field(ax1,ax2,fc_params,field)
      
   else:
      raise Exception('Ouch! Instrument unknown ...')   

# ----------------------------------------------------------------------------------------
def get_total_right_radius(fc_params, right_radius):
   '''
   A function to compute the total right radius given the offsets within the OB.
   
   Args:
      fc_params: dictionnary of OB parameters
      right_radius: default right radius
   Return: 
      right_radius: updated right_radius
   '''
   
   # loop through all the Science template, and find the largest offset
   for n in range(fc_params['n_sci']):
      off1 = fc_params['sci%i'%(n+1)]['off1']
      off2 = fc_params['sci%i'%(n+1)]['off1']
      offs_abs = [ np.array([np.sum(off1[:i]) for i in range(len(off1))]) + fc_params['acq']['bos_ra'], 
                   np.array([np.sum(off2[:i]) for i in range(len(off2))]) + fc_params['acq']['bos_dec'],
                  ]
      # Radius of all steps
      rad_obs = np.sqrt( (offs_abs[0])**2 + (offs_abs[1])**2)
      # Select the largets one
      right_radius = np.max([right_radius,np.max(rad_obs)+fcm_xshooter.right_radius])
   
   return right_radius

# ----------------------------------------------------------------------------------------
def get_chart_radius(fc_params):
   '''
   Defines the radius of the right and left plots for different instruments.
   
   Args:
      fc_params: the parameters of the finding charts (a dictionary)
   Returns:
      (left_radius, right_radius): a tuple with the chart radius in arcsec.
   '''
   
   if fc_params['inst'] == 'MUSE':
   
      # For the right hand-side, find the largest possible offset, to make sure even far-away
      # sky fields are visible
      if fc_params['n_sci'] == 0:
         right_radius = fcm_muse.right_radius
      else:
         
         right_radius = copy.deepcopy(fcm_muse.right_radius)
         # Account for the OB offsets
         right_radius = get_total_right_radius(fc_params, right_radius)
   
      left_radius = fcm_muse.left_radius(fc_params['ins_mode']) + 0.5* np.sqrt( (fc_params['acq']['bos_ra']/2.)**2 + 
                                     (fc_params['acq']['bos_dec']/2.)**2)
      
   elif fc_params['inst'] == 'HAWKI':
      
      # For the right-hand side, adjust the radius to fit in the largest offsests
      if fc_params['n_sci'] == 0:
         right_radius = fcm_hawki.right_radius
      else:
         right_radius = copy.deepcopy(fcm_hawki.right_radius)
         # Account for the OB offsets
         right_radius = get_total_right_radius(fc_params, right_radius)
      
      # Keep Left hand-side focused on acq field
      left_radius = fcm_hawki.left_radius
   
   elif fc_params['inst'] == 'XSHOOTER':
   
      # For the right-hand side, adjust the radius to fit in the largest offsests
      if fc_params['n_sci'] == 0:
         right_radius = fcm_xshooter.right_radius
      else:
         right_radius = copy.deepcopy(fcm_xshooter.right_radius)
         # Account for the OB offsets
         right_radius = get_total_right_radius(fc_params, right_radius)

      # Keep Left hand-side focused on acq field
      left_radius = fcm_xshooter.left_radius + 0.5* np.sqrt( (fc_params['acq']['bos_ra']/2.)**2 + 
                                     (fc_params['acq']['bos_dec']/2.)**2)
   
   else:
      raise Exception('Ouch! Instrument unknown ...')
   
   return (left_radius, right_radius)

# ----------------------------------------------------------------------------------------
def get_GS_outer_radius(fc_params):
   '''
   Returns the outer search radius for different instruments.
   
   Args:
      fc_params: the parameters of the finding charts (a dictionary)
   Returns:
      radius: the radius in arcsec.
   '''
   
   if fc_params['inst'] in ['MUSE', 'HAWKI']:
      return fcm_m.outer_GS_Nas
   elif fc_params['inst'] in ['XSHOOTER']:
      return fcm_m.outer_GS_Cas
   else:
      raise Exception('Ouch! Instrument unknown ...')

# ----------------------------------------------------------------------------------------
def get_gaia_image_params(fc_params):
   '''
   Returns the outer search radius for different instruments.
   
   Args:
      fc_params: the parameters of the finding charts (a dictionary)
   Returns:
      (sampling,fwhm): the sampling and fwhm to build a fake Gaia image 
   '''
   
   if fc_params['inst'] == 'MUSE':
      if 'NFM' in fc_params['ins_mode']:
         return (0.025*u.arcsec, 0.08*u.arcsec)
      else: 
         return (0.3*u.arcsec, 0.6*u.arcsec)
         
   elif fc_params['inst'] == 'HAWKI':
      return (0.3*u.arcsec, 1.0*u.arcsec) 
      
   elif fc_params['inst'] == 'XSHOOTER':
      return (0.3*u.arcsec, 1.0*u.arcsec)
      
   else:
      raise Exception('Ouch! Instrument unknown ...')

# ----------------------------------------------------------------------------------------
def get_chart_center(fc_params):
   '''
   Defines the center of the right and left plots for different instruments.
   
   Args:
      fc_params: the parameters of the finding charts (a dictionary)
   Returns:
      (left_center, right_center): a tuple with the chart center as a SkyCoord value.
   '''
   
   if fc_params['inst'] == 'MUSE':
      # Place the field center between the target and the acquisiton
      center = fcm_t.offset_coord(fc_params['target'], 
                                  delta_ra=fc_params['acq']['bos_ra']/2.*u.arcsec, 
                                  delta_dec=fc_params['acq']['bos_dec']/2.*u.arcsec)
   
      return (center, center)
      
   elif fc_params['inst'] == 'HAWKI':
      # Place the chart center on the acquisition field
      center = fcm_t.offset_coord(fc_params['target'], 
                                  delta_ra= fc_params['acq']['bos_ra']*u.arcsec, 
                                  delta_dec= fc_params['acq']['bos_dec']*u.arcsec)
   
      return (center, center) 
   
   elif fc_params['inst'] == 'XSHOOTER':
      # Place the field center between the target and the acquisiton
      # Warning, remember that with XSHOOTER, the offset are from the BOS to the target !!!
      center = fcm_t.offset_coord(fc_params['target'], 
                                  delta_ra=-fc_params['acq']['bos_ra']/2.*u.arcsec, 
                                  delta_dec=-fc_params['acq']['bos_dec']/2.*u.arcsec)  
                                  
      return (center, center)                             
   
   else:
      raise Exception('Ouch! Instrument unknown ...')

# ----------------------------------------------------------------------------------------
def get_scalebar(inst, ins_mode = None):
   '''
   Sets the scalebar of the chart, as a function of the instrument
   
   Args:
      inst (string): the instrument name
   Returns:
      (scl, scll): tuple with the scalebar length (in degrees) and label
   '''
   
   if inst == 'MUSE':
      if ins_mode is None:
         raise Exception('Ouch! For MUSE, I need to know the mode !')
      elif 'NFM' in ins_mode:
         return (1./3600, '1$^{\prime\prime}$')
      else:
         return (20./3600, '20$^{\prime\prime}$')
   elif inst == 'HAWKI':
      return (60./3600, '1$^{\prime}$')
   elif inst == 'XSHOOTER':
      return (20./3600, '20$^{\prime\prime}$')
   else:
      raise Exception('Ouch! Instrument unknown ...')

# ----------------------------------------------------------------------------------------
def get_pmin(survey):
   '''
   A function to get the pmin value for the aplpy chart, as a function of the instrument.
   
   Args:
      fc_params: the parameters of the finding charts (a dictionary)
   Returns:
      pmin (float): the pmin value
   '''
   
   if survey in ['2MASS-J', '2MASS-H', '2MASS-K']:
      return 30.
   
   else:
      return 10.

# ----------------------------------------------------------------------------------------
def get_inner_GS_search(fc_params):
   '''
   A function to get theinner search radius for Gudie Stars, as a function of the instrument.
   
   Args:
      fc_params: the parameters of the finding charts (a dictionary)
   Returns:
      radius (float): the radius value, in arcsec
   '''
   
   if fc_params['inst'] == 'MUSE':
      return fcm_muse.inner_GS_search
      
   elif fc_params['inst'] == 'HAWKI':
      return fcm_hawki.inner_GS_search
      
   elif fc_params['inst'] == 'XSHOOTER':
      return fcm_xshooter.inner_GS_search
      
   else:
      raise Exception('Ouch! Instrument unknown ...')
# ----------------------------------------------------------------------------------------
def get_target_from_ephem(fc_params, ephemeris):
   '''
   A function to extract the target coordinates and other useful stuff from a PAF
   ephemeris file.
    
   Args:
      fc_params: the parameters of the finding charts (a dictionary)
      ephemeris: the content of a PAF files as a list of strings from readlines()
   Returns:
      fc_params: the updated parameters of the finding charts (a dictionary)
    
    '''
   # First, extract all the ephemeris points, store them in a list
   # Time, Ra, Dec, pmra, pmdec
   # This VERY MUCH ASSUMES a PAF file !!!
   ephem = [list( line.split('"')[1].split(',')[i] for i in [0,2,3,4,5]) 
            for line in ephemeris if line[:16]=='INS.EPHEM.RECORD'] 
      
   # Extract all the times, convert them to datetime values. UTC = PAF default!
   times = [ dup.parse(point[0]+' UTC') for point in ephem]
      
   # Here, let's do some sanity check:
   if (fcm_m.obsdate>times[-1]) or (fcm_m.obsdate<times[0]):
      raise Exception('Ouch! Specified obstime outside of Ephemeris range: %s - %s' % 
                       (times[0].strftime('%Y-%m-%d %H:%M:%S %Z'),
                        times[-1].strftime('%Y-%m-%d %H:%M:%S %Z'))) 
      
   # Find the index of the closest ephemeris point
   closest = np.argmin(np.abs(np.array(times)-fcm_m.obsdate))
      
   # Check if I will be accurate (or not)
   if np.abs(times[closest]-fcm_m.obsdate) > timedelta(minutes=30):
      warnings.warn('Observing time %.1f away from closest ephemeris point' % 
                     (np.abs(times[closest]-fcm_m.obsdate).total_seconds()/60. ))
      
   closest_target = SkyCoord(ra=ephem[closest][1], dec=ephem[closest][2], 
                            unit=(u.hourangle, u.degree),
                            pm_ra_cosdec = float(ephem[closest][3]) *u.arcsec/u.second,
                            pm_dec = float(ephem[closest][3])*u.arcsec/u.second, 
                            obstime= Time(times[closest]),
                            distance=fcm_m.ephem_d,
                            )
   # Derive the target by propagating proper motions from 
   fc_params['target'] = closest_target.apply_space_motion(new_obstime = Time(fcm_m.obsdate))
      
   # Now, I also want to show all the points in the ephemeris file with a fcm_m.ephem_range time interval
   # Find the index and associated times of all the near-by ephemeris entries
   # The advantage, for these, is that I make no errors by assuming a distance

   # For plotting purposes, deal with past and future ephemeris entry points separately   
   close_past = (timedelta(seconds=0) <= -(np.array(times)-fcm_m.obsdate)) * (-(np.array(times)-fcm_m.obsdate) <= timedelta(seconds = fcm_m.ephem_range.to(u.second).value))
      
   close_future = (timedelta(seconds=0) <= (np.array(times)-fcm_m.obsdate)) * ((np.array(times)-fcm_m.obsdate) <= timedelta(seconds = fcm_m.ephem_range.to(u.second).value))
    
   # Now store the points as proper SkyCoord entries  
   fc_params['ephem_points_past'] = [SkyCoord(ra=item[1], dec=item[2], 
                                        unit=(u.hourangle, u.degree),
                                        pm_ra_cosdec = float(item[3]) *u.arcsec/u.second,
                                        pm_dec = float(item[3])*u.arcsec/u.second, 
                                        obstime= Time(dup.parse(item[0]+' UTC')),
                                        distance=fcm_m.ephem_d,
                                        )  for item in list(np.array(ephem)[close_past])]
   past_times = np.array(times)[close_past] 
   # Also store the time delta, for legend purposes.                                 
   fc_params['ephem_past_delta'] = np.median(past_times[1:]-past_times[:-1]) 
       
   # Idem for the future points.                                  
   fc_params['ephem_points_future'] = [SkyCoord(ra=item[1], dec=item[2], 
                                        unit=(u.hourangle, u.degree),
                                        pm_ra_cosdec = float(item[3]) *u.arcsec/u.second,
                                        pm_dec = float(item[3])*u.arcsec/u.second, 
                                        obstime= Time(dup.parse(item[0]+' UTC')),
                                        distance=fcm_m.ephem_d,
                                        )  for item in list(np.array(ephem)[close_future])]
      
   future_times = np.array(times)[close_future]    
   # For legend purposes.                              
   fc_params['ephem_future_delta'] = np.median(future_times[1:]-future_times[:-1]) 
   
   return fc_params

      
# ----------------------------------------------------------------------------------------
def get_p2fcdata(obID, api):
   ''' 
   Extract the required info for creating a finding chart from the ob and templates
   obtained via p2api.
   
   Args:
      obID: the OB Id, as an integer.
      api: a p2api.ApiConnection() object 
      
   Returns:
      a dictionnary containing the OB parameters
   '''
   
   # Fetch the OB
   ob, obVersion = api.getOB(obID)
   
   # Check for an ephemeris file ...
   # This will download a file, even if there is no ephemeris file upload by the user
   ephem_fn = os.path.join(fcm_m.data_loc,'ephem_%i.txt'%(obID))
   api.getEphemerisFile(obID, ephem_fn)
   
   # I have no choice but to open it, check what's inside ...
   f = open(ephem_fn)
   ephemeris = f.readlines()
   f.close()
   
   if len(ephemeris) == 0:
   # Ok, there's no ephemeris file attached to the OB ... clean-up my mess
      os.remove(ephem_fn)
      ephem_fn = None
   else:
      print('   I found an ephemeris file, and stored it: %s' % ephem_fn)
   
   # Start filling my own dictionary of useful finding chart parameters
   fc_params = {}
   
   fc_params['inst'] = ob['instrument']
   fc_params['ins_mode'] = None # Create a general entry, to be tweaked later on if needed
   fc_params['ob_name'] = ob['name']
   
   fc_params['ob_id'] = obID
   
   # Get the PI name and progId
   runId = ob['runId']
   run, _ = api.getRun(runId)
   fc_params['pi'] = run['pi']['lastName']
   fc_params['prog_id'] = run['progId']
   
   # Also create a list where I keep track of anything special about the OB
   # possible flag include 'moving_target', 'parallactic_angle'
   fc_params['tags'] = []
   
   # Need to deal with the target coordinates. Ephemeris file, or not ?
   if ephem_fn is None:
    
      # Issue #6: https://github.com/fpavogt/fcmaker/issues/6
      # Make sure the equinox comes with a "J"
      if not(ob['target']['equinox'][0] == 'J'):
         warnings.warn('Equinox %s invalid: it should start with a "J". Using %s instead.' 
                       % (ob['target']['equinox'], 'J'+ob['target']['equinox']) )
         ob['target']['equinox'] = 'J'+ob['target']['equinox']
    
      # Ok, just a normal target ... extract the required info
      tc = SkyCoord(ra = ob['target']['ra'],
                    dec = ob['target']['dec'], 
                    obstime = Time(ob['target']['epoch'],format='decimalyear'),
                    equinox = ob['target']['equinox'], 
                    frame = 'icrs',unit=(u.hourangle, u.deg), 
                    pm_ra_cosdec = ob['target']['properMotionRa']*u.arcsec/u.yr,
                    pm_dec = ob['target']['properMotionDec']*u.arcsec/u.yr,
                    # I must specify a generic distance to the target,
                    # if I want to later on propagate the proper motions
                    distance=100*u.pc,  
                    )
                                  
      # Propagate the proper motion using astropy v3.0
      fc_params['target'] = tc.apply_space_motion(new_obstime = Time(fcm_m.obsdate))  
      
      fc_params['ephem_points_past'] = []                               
      fc_params['ephem_points_future'] = []
      # Flag it if this is a moving target 
      if np.sqrt(ob['target']['properMotionRa']**2 + ob['target']['properMotionDec']**2) > 0.0:
         fc_params['tags'] += ['moving_target']

   else:
   
      # Now deal with this ephemeris file
      fc_params = get_target_from_ephem(fc_params, ephemeris)
      fc_params['tags'] += ['moving_target']
   
     
   # Extract the acquisition and observations parameters for the support instruments
   if ob['instrument'] == 'MUSE':
      return fcm_muse.get_p2fcdata_muse(fc_params, ob, api)
      
   elif ob['instrument'] == 'HAWKI':
      return fcm_hawki.get_p2fcdata_hawki(fc_params, ob, api)
   
   elif ob['instrument'] == 'XSHOOTER':
      return fcm_xshooter.get_p2fcdata_xshooter(fc_params, ob, api)
      
   else:
      raise Exception('%s finding charts not (yet?) supported.' % (inst))
   
# ----------------------------------------------------------------------------------------
def get_localfcdata(inpars):
   ''' 
   Extract the required info for creating a finding chart from the ob and templates
   obtained from a local parameter file.
   
   Args:
      inpars: a dictionnary.
      
   Returns:
      a dictionnary containing the OB parameters
   '''

   # Deal with the ephemeris ...
   fc_params = {}

   #fc_params = {}
   fc_params['ob_name'] = inpars['ob_name']
   
   if type(inpars['ob_id']) == int:
      fc_params['ob_id'] = inpars['ob_id']
   else:
      fc_params['ob_id'] = -1
      
   fc_params['pi'] = inpars['pi']
   fc_params['prog_id'] = inpars['prog_id']
   fc_params['inst'] = inpars['inst']
   fc_params['ins_mode'] = None # Create a general entry, to be tweaked later on if needed
   
   # Also create a list where I keep track of anything special about the OB
   # possible flag include 'moving_target', 'parallactic_angle'
   fc_params['tags'] = []

   # Some other day ...
   #fc_params['ephem_points_past'] = []                               
   #fc_params['ephem_points_future'] = [] 
   
   if inpars['ephemeris'] is None:
      # Ok, just a normal target ... extract the required info
      tc = SkyCoord( ra = inpars['ra'],
                     dec = inpars['dec'], 
                     obstime = Time(inpars['epoch'],format='decimalyear'),
                     equinox = inpars['equinox'], 
                     frame = 'icrs',unit=(u.hourangle, u.deg), 
                     pm_ra_cosdec = inpars['pmra']*u.arcsec/u.yr,
                     pm_dec = inpars['pmdec']*u.arcsec/u.yr,
                     # I must specify a generic distance to the target,
                     # if I want to later on propagate the proper motions
                     distance=fcm_m.default_pm_d,  
                     )
                                  
      # Propagate the proper motion using astropy v3.0
      fc_params['target'] = tc.apply_space_motion(new_obstime = Time(fcm_m.obsdate))  
      # Flag it as a moving target if needed 
      if np.sqrt(inpars['pmra']**2 + inpars['pmdec']**2) > 0.0:
         fc_params['tags'] += ['moving_target'] 
                                   
      fc_params['ephem_points_past'] = []                               
      fc_params['ephem_points_future'] = []   
      
   else:
      ephem_fn = os.path.join('.',inpars['ephemeris'])
      
      if not(os.path.isfile(ephem_fn)):
         raise Exception('Ouch! No ephemeris file found at: %s'% (ephem_fn))
      
      # Ope the ephemeris file
      f = open(ephem_fn)
      ephemeris = f.readlines()
      f.close()
   
      # Now deal with this ephemeris file
      fc_params = get_target_from_ephem(fc_params, ephemeris)
      fc_params['tags'] += ['moving_target']
      
   # Extract the acquisition and observations parameters for the support instruments
   if inpars['inst'] == 'MUSE':
      return fcm_muse.get_localfcdata_muse(fc_params, inpars)
   elif inpars['inst'] == 'HAWKI':
      return fcm_hawki.get_localfcdata_hawki(fc_params, inpars) 
   elif inpars['inst'] == 'XSHOOTER':
      return fcm_xshooter.get_localfcdata_xshooter(fc_params, inpars)
   else:
      raise Exception('%s finding charts not (yet?) supported.' % (inpars['inst']))
     
# ------------------------------------------------------------------------------
      