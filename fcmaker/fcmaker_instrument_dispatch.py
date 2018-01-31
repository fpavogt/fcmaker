# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------------------

import numpy as np
import warnings

from astropy.coordinates.sky_coordinate import SkyCoord
from astropy import units as u
import copy

from . import fcmaker_metadata as fcm_m
from . import fcmaker_tools as fcm_t
from . import fcmaker_muse as fcm_muse
from . import fcmaker_hawki as fcm_hawki


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
      return fcm_muse.bk_image
      
   elif fc_params['inst'] == 'HAWKI':
      # For HAWKI, in case of J, H or K, then use the corresponding 2MASS image
      if fc_params['acq']['filter'] in ['J','H','Ks']:
         return '2MASS-'+fc_params['acq']['filter'][0]
      else:
         return fcm_hawki.bk_image

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
      
      
      #return (ax1, ax2)
      
   else:
      raise Exception('Ouch! Instrument unknown ...')   

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
         # loop through all the Science template, and find the largest offset
         for n in range(fc_params['n_sci']):
            off1 = fc_params['sci%i'%(n+1)]['off1']
            off2 = fc_params['sci%i'%(n+1)]['off1']
            offs_abs = [ np.array([np.sum(off1[:i]) for i in range(len(off1))]), 
                         np.array([np.sum(off2[:i]) for i in range(len(off2))]),
                       ]
            # Radius of all steps
            rad_obs = np.sqrt( (offs_abs[0])**2 + (offs_abs[1])**2)
            # Select the largets one
            right_radius = np.max([right_radius,np.max(rad_obs)+60.])
         
      # Ok, I want to have a "zoom-in version of the plot on the left. 
      # In case of big offsets, this can be a problem. So let's make this a bit larger
      # in those cases
      left_radius = np.max([fcm_muse.left_radius,
                            np.sqrt( (fc_params['acq']['bos_ra']/2.)**2 + 
                                     (fc_params['acq']['bos_dec']/2.)**2)+35.])
   
   elif fc_params['inst'] == 'HAWKI':
      
      # For the right-hand side, adjust the radius to fit in the largest offsests
      if fc_params['n_sci'] == 0:
         right_radius = fcm_hawki.right_radius
      else:
         right_radius = copy.deepcopy(fcm_hawki.right_radius)
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
            right_radius = np.max([right_radius,np.max(rad_obs)+240.])
      
      # Keep Left hand-side focused on acq field
      left_radius = fcm_hawki.left_radius
   
   else:
      raise Exception('Ouch! Instrument unknown ...')
   
   return (left_radius, right_radius)

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
   
   else:
      raise Exception('Ouch! Instrument unknown ...')

# ----------------------------------------------------------------------------------------
def get_scalebar(inst):
   '''
   Sets the scalebar of the chart, as a function of the instrument
   
   Args:
      inst (string): the instrument name
   Returns:
      (scl, scll): tuple with the scalebar length (in degrees) and label
   '''
   
   if inst == 'MUSE':
      return (30./3600, '30$^{\prime\prime}$')
   elif inst == 'HAWKI':
      return (60./3600, '1$^{\prime}$')
      
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
      
   else:
      raise Exception('Ouch! Instrument unknown ...')
      
# ----------------------------------------------------------------------------------------
def get_p2fcdata(obID,api):
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
   inst = ob['instrument']
    
   # Extract the acquisition and observations parameters for the support instruments
   if inst == 'MUSE':
      return fcm_muse.get_p2fcdata_muse(ob, api)
      
   elif inst == 'HAWKI':
      return fcm_hawki.get_p2fcdata_hawki(ob, api)
      
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

   # Extract the acquisition and observations parameters for the support instruments
   if inpars['inst'] == 'MUSE':
      return fcm_muse.get_localfcdata_muse(inpars)
   elif inpars['inst'] == 'HAWKI':
      return fcm_hawki.get_localfcdata_hawki(inpars)  
   else:
      raise Exception('%s finding charts not (yet?) supported.' % (inst))
     
# ------------------------------------------------------------------------------
      