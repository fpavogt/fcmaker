# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------------------

import numpy as np
from astropy.coordinates.sky_coordinate import SkyCoord
import astropy.units as u
from astropy.time import Time
from astroplan import Observer
import warnings
import copy

import matplotlib.lines as mlines

from . import fcmaker_tools as fcm_t
from . import fcmaker_metadata as fcm_m

'''
fcmaker: a Python module to automatically create finding charts for ESO OBs in p2.\n
Copyright (C) 2017-2019,  F.P.A. Vogt
--- oOo ---
This file contains tools related to the ESPRESSO instrument.
Created September 2018, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''


# ---------------- Define the ESPRESSO field of views -------------------------------------

espresso_acq_fov = 17*u.arcsec # in radius

espresso_fiber_sep = 7.5*u.arcsec # Distance between fiber A and fiber B

espresso_fov = {'HR': 0.5*u.arcsec, # in radius ! Manual p11
                'UHR': 0.25*u.arcsec,
                'MR': 0.5*u.arcsec,
               }

# ----------------------------------------------------------------------------------------

# The default background image from SkyView
bk_image = 'Gaia'

# The radius of the charts
right_radius = 720. # in arcsec

left_radius = 30. # in arcsec
      
inner_GS_search = 120. # inner limit to find Guide Stars # TO BE VALIDATED!

# List the supported ESPRESSO observing templates
espresso_templates = [# HR
                      'ESPRESSO_singleHR_acq_obj',
                      'ESPRESSO_singleHR_obs',
                      # UHR
                      'ESPRESSO_singleUHR_acq_obj', 
                      'ESPRESSO_singleUHR_obs',
                      #multiMR
                      'ESPRESSO_multiMR_acq_obj',
                      'ESPRESSO_multiMR_obs',
                      ]

# The acquisition and science parameters that matter for the finding charts.

espresso_acq_params = {'bos_ra': 0,      # blind offset RA
                       'bos_dec': 0,     # blind offset Dec
                       'is_gs': False,   # Guide star defined by user ?
                       'gs': None,       # gs coord as SkyCoord
                  }    
         
espresso_sci_params = {'ins_mode':None, # Can be 'HR', 'UHR', 'MR'
                       'calsource_B': 'DARK', # Can be 'DARK','FPCS','SKY', 'THAR2'
                       }

# ----------------------------------------------------------------------------------------
def get_p2fcdata_espresso(fc_params, ob, api):
   '''
   Extracts all the important info to build a finding chart from a given ESPRESSO OB from p2.
   
   Args:
      fc_params: dictionnary of finding chart parameters
      ob: an api.getOB() object
      api: a p2api.ApiConnection() object 
   Returns:
       A dictionnary containing the OB parameters
   '''
   
   obId = ob['obId']
   
   # Fetch all the associated templates
   templates, templatesVersion = api.getTemplates(obId)
   
   fc_params['n_sci'] = 0 # The number of Science template in the OB 
   
   # First of all, check if I have an acquisition template - required to set the 
   # instrument mode, etc ...
   is_acq = False
   for t in templates:
      if t['type'] == 'acquisition':
         is_acq = True
   
   if  not(is_acq):
      raise Exception('Ouch! I found no acquisition template (required)!')
   
   # Ok, from now on, I can assume that the OB does contain an acq template.
   
   # First, let us loop through all the templates in the OB, to see which ones are supported
   for (n,t) in enumerate(templates):
      # If this is a weird template, issue a warning, and continue
      if not(t['templateName'] in espresso_templates):
         warnings.warn('Template %s is not supported by fcmaker. Skipping it.' % (t['templateName']))
         continue

      # Ok, this is a template I know about. Let's check which it is.
      # p2api returns the templates in the order defined inside p2 by the user. 
      # I.e. it places the acquisition first, even if it was not create first.
      if t['type'] == 'acquisition':
         
         # There should only ever be one acquisition ... let's make sure of this ...
         try:
            fc_params['acq'] # Not yet created, should create an error
         except:
            pass # This is the expected route
         else: # fc_params['acq'] already exists - this should not be! 
            raise Exception('Ouch! There can only be one acquisition template per OB!')
            
         # Fill the acq parameters
         fc_params['acq'] = copy.deepcopy(espresso_acq_params)
         
         if 'UHR' in t['templateName']:
            fc_params['acq']['ins_mode'] = 'UHR'
         elif 'HR' in t['templateName']:
            fc_params['acq']['ins_mode'] = 'HR'
         elif 'MR' in t['templateName']:
            fc_params['acq']['ins_mode'] = 'MR' 
         else:
            raise Exception('Ouch! This error is impossible!')
         
         
         tpl,tplVersion = api.getTemplate(obId, t['templateId'])
         for param in tpl['parameters']:
                      
            # Guide Star
            if param['name'] == 'TEL.AG.GUIDESTAR':
               if param['value'] in ['NONE','CATALOGUE']:              
                  fc_params['acq']['is_gs'] = False
               else:
                  fc_params['acq']['is_gs'] = True
                             
            if param['name'] == 'TEL.GS1.ALPHA':             
               gs_ra = param['value']       
            if param['name'] == 'TEL.GS1.DELTA':             
               gs_dec = param['value']
             
            # Blind offset           
            if param['name'] == 'TEL.TARG.OFFSETALPHA':             
               fc_params['acq']['bos_ra'] = param['value']  
            if param['name'] == 'TEL.TARG.OFFSETDELTA':             
               fc_params['acq']['bos_dec'] = param['value']    
            
         # Store the GS and TTS as SkyCoords
         if fc_params['acq']['is_gs']:
            fc_params['acq']['gs'] = SkyCoord(gs_ra, gs_dec, 
                                       obstime = fcm_m.obsdate, 
                                       #equinox=ob['target']['equinox'], 
                                       frame='icrs',unit=(u.hourangle, u.deg))
      
      
      elif t['type'] in ['science','calib']:
         
         # Ok, I found one more Science template
         fc_params['n_sci'] += 1
         
         temp_name = 'sci%i' %(fc_params['n_sci'])
         
         # Store all the relevant parameters in a dictionary
         fc_params[temp_name] = copy.deepcopy(espresso_sci_params)
         
         tpl,tplVersion = api.getTemplate(obId, t['templateId'])
         
         if 'UHR' in t['templateName']:
            fc_params[temp_name]['ins_mode'] = 'UHR'
         elif 'HR' in t['templateName']:
            fc_params[temp_name]['ins_mode'] = 'HR'
         elif 'MR' in t['templateName']:
            fc_params[temp_name]['ins_name'] = 'MR'
         else:
            raise Exception('Ouch! This error is impossible!')
            
         for param in tpl['parameters']:

            if param['name'] == 'SEQ.CALSOURCEB':             
               fc_params[temp_name]['calsource_B'] = param['value']
                      
   return fc_params

# ------------------------------------------------------------------------------
def get_localfcdata_espresso(fc_params,inpars):
   '''
   Extracts all the important info to build a finding chart from a given ESPRESSO OB defined
   locally.
   
   Args:
      inpars: A dictionnary containing the OB parameters
   Returns:
       A dictionnary containing the ESPRESSO OB parameters
   '''                                    
         
   # Acquisition
   fc_params['acq'] = copy.deepcopy(espresso_acq_params)
   fc_params['acq']['ins_mode'] = inpars['ins_mode']
   fc_params['acq']['bos_ra'] = inpars['bos_ra']
   fc_params['acq']['bos_dec'] = inpars['bos_dec']
   fc_params['acq']['is_gs'] = inpars['is_gs']
   fc_params['acq']['gs'] = SkyCoord(inpars['gs_ra'],inpars['gs_dec'], 
                                     frame = 'icrs', 
                                     obstime = Time(fcm_m.obsdate),
                                     equinox = 'J2000')

   # "Observation" ... for ESPRESSO, just to known if Fiber B is on SKY or not
   fc_params['n_sci'] = 1
   fc_params['sci1'] = copy.deepcopy(espresso_sci_params)
   fc_params['sci1']['ins_mode'] = inpars['ins_mode']
   fc_params['sci1']['calsource_B'] = inpars['calsource_B']

   return fc_params
   
# ------------------------------------------------------------------------------
def get_fields_dict(fc_params):
   '''
   Create a dictionnary with the basic info required to draw the fields, given 
   certain OB parameters.
   
   Args:
      fc_params: A dictionnary containing the ESPRESSO OB parameters
   Returns:
       A dictionary containing the plotting parameters for each field.
   
   '''
   
   # Some values to keep track of the ra, dec and pa offsets over the course of the OB.
   delta_ra = 0
   delta_dec = 0

   # Define the center of the finding chart to be on the acquisition field.
   target = fc_params['target']  
   
   # Observing fields : build a dictionary for all of them
   fields = {}

   # First, the acquisiton frame 
   counter = 1
   acq_field =  fcm_t.offset_coord(target, 
                                   delta_ra = fc_params['acq']['bos_ra']*u.arcsec, 
                                   delta_dec = fc_params['acq']['bos_dec']*u.arcsec)
   
   # Start filling the dictionary of fields
   fields[1] = [fc_params['inst'],
                fc_params['acq']['ins_mode'],
                 acq_field, # Field central coordinates
                'Acq', # Nature of the field ('Acq', 'Target')
                espresso_acq_fov, # The size of the FoV 
                ]
      
   # Include the Target as a "field"
   counter+=1
   fields[counter] = [fc_params['inst'],
                      fc_params['acq']['ins_mode'],
                      target, # Field central coordinates
                      'Target', # Nature of the field ('Acq', 'Target')
                      espresso_acq_fov, # The ins_mode
                      ]  
   
   
   # Now, loop through all the Science templates I have
   # I do this only to check if anyone of them has the Fiber B on sky.
   # All coordinates remain exactly the same as the Target.
   for n in range(fc_params['n_sci']):
      
      # What is the name of this specific template ?
      temp_name = 'sci%i' %(n+1)
      
      # Ok, what is the calsource B ?
      calsource_B = fc_params[temp_name]['calsource_B']
      
      if calsource_B == 'SKY':
         
         # Only store the "field" if I am on -sky. Else, just keep the target.
         counter+=1
         fields[counter] = [fc_params['inst'],
                            fc_params['acq']['ins_mode'],
                            target, # Field central coordinates
                            'calsource', # Nature of the field ('Acq', 'Target')
                            espresso_acq_fov, # The ins_mode
                           ] 
   
   return fields

# ----------------------------------------------------------------------------------------
def plot_field(ax1, ax2, fc_params, field):
   '''
   The specific ESPRESSO function that draws a specific observation field.
   
   Args:
      ax1,ax2: the left and right plots 'axis'
      fc_params: the dictionnary of parameters for the OB
      field: the specific field list of parameters
   '''
   
   skins = {'Acq':{'c':'darkorchid', 'lwm':2, 'zorder':10, 'marker':fcm_t.crosshair(pa=45), 
                   'lw':1.5, 'ms':500, 'ls':'-', 'mc':'darkorchid'},
            'Target': {'c':'darkorange', 'lwm':1, 'zorder':5, 'marker':'o', 
                       'lw':1.0, 'ms':250, 'ls':'--', 'mc':'darkorange'},
            'calsource': {'c':'darkorange', 'lwm':1, 'zorder':5, 'marker':'None', 
                          'lw':1.0, 'ms':250, 'ls':':', 'mc':'darkorange'},
            }
   
   this_coords = [field[2].ra, field[2].dec]
   
   # For the acquisition, mark the center with a crosshair
   if field[3] == 'Acq':
      ax1.show_markers(this_coords[0].deg, this_coords[1].deg, 
                       marker=skins[field[3]]['marker'],
                       edgecolor=skins[field[3]]['mc'],
                       s=skins[field[3]]['ms'], 
                       linewidth=skins[field[3]]['lwm'], 
                       zorder=skins[field[3]]['zorder'],
                       )
                         
   # For the target, show the actual fiber size
   if field[3] == 'Target':
      ax1.show_circles([this_coords[0].deg],
                       [this_coords[1].deg],
                       [espresso_fov[field[1]].to(u.degree).value],
                       color=skins[field[3]]['mc'], 
                       lw=skins[field[3]]['lw'], 
                       zorder=skins[field[3]]['zorder'])
   
   # Then show the full acq camera field-of-view   
   for ax in [ax1,ax2]:
       
      if field[3] in ['Acq','Target']:
  
         # Also show the acquisition camera field-of-view
         ax.show_circles([this_coords[0].deg],
                         [this_coords[1].deg],
                         [espresso_acq_fov.to(u.degree).value],
                         color=skins[field[3]]['mc'], 
                         lw=skins[field[3]]['lw'], 
                         linestyle=skins[field[3]]['ls'],
                         zorder=skins[field[3]]['zorder']) 
                             
         # and the GS validity area
         ax.show_circles([this_coords[0].deg],
                      [this_coords[1].deg],
                      [((fcm_m.outer_GS_Nas*u.arcsec).to(u.degree)).value,],
                       color='k', lw=0.5, 
                       zorder=skins[field[3]]['zorder']
                     )
      
      # In case I have Fiber B on-sky, show it's possible location                          
      elif field[3] in ['calsource']:

         ax.show_circles([this_coords[0].deg],
                         [this_coords[1].deg],
                         [espresso_fiber_sep.to(u.degree).value],
                          color=skins[field[3]]['mc'], 
                          linestyle=skins[field[3]]['ls'],
                          lw=skins[field[3]]['lw'], 
                          zorder=skins[field[3]]['zorder'])
                                    
                   
      # Show the Guide Star (if there is an acq template present)
      if (fc_params['acq']['is_gs']) and (field[3] == 'Acq'):
         ax.show_markers(fc_params['acq']['gs'].ra,fc_params['acq']['gs'].dec, 
                         marker=fcm_t.crosshair(pa=45),edgecolor='crimson',
                         s=500, linewidth=1.5)
      
         if fcm_m.fcm_usetex:
            lab = r'\textbf{GS}'
         else: 
            lab = r'GS'
         
         if ax == ax2:                  
            # Only add the name to the large plot
            ax2.add_label(fc_params['acq']['gs'].ra.deg,fc_params['acq']['gs'].dec.deg+70./3600,lab, 
                    verticalalignment='center', 
                    horizontalalignment='center',size=12,color='k',
                    bbox=dict(facecolor='w',ec='k', alpha=0.6))    
      
            # Check if G is compatible with GS area for this instrument
            if ((field[2].separation(fc_params['acq']['gs']) > (fcm_m.outer_GS_Nas*u.arcsec)) or \
                (field[2].separation(fc_params['acq']['gs']) < (inner_GS_search*u.arcsec))):
                
               if fcm_m.fcm_usetex:
                  lab = r'\textbf{!}'
               else:
                  lab = '!'
                  
               ax.add_label(fc_params['acq']['gs'].ra.deg,fc_params['acq']['gs'].dec.deg,lab, 
                            verticalalignment='center', 
                            horizontalalignment='center',size=12,color='crimson',
                            bbox=dict(boxstyle="circle,pad=0.17", facecolor='w',ec='crimson', alpha=1))
          
   
      # Finally, also add a legend for clarity
      acq_legend = mlines.Line2D([], [], 
                                  color=skins['Acq']['c'],
                                  markerfacecolor='None',
                                  markeredgecolor=skins['Acq']['mc'], 
                                  marker=skins['Acq']['marker'],
                                  linestyle=skins['Acq']['ls'],
                                  linewidth=skins['Acq']['lw'],
                                  markersize=10, label='Acq.')
      target_legend = mlines.Line2D([], [],
                                    color=skins['Target']['c'],
                                    markerfacecolor='None', 
                                    markeredgecolor=skins['Target']['mc'], 
                                    linestyle=skins['Target']['ls'],
                                    linewidth=skins['Target']['lw'],
                                    marker=skins['Target']['marker'],
                                    markersize=10, label='Target')
      calsource_legend = mlines.Line2D([], [], 
                               color=skins['calsource']['c'],
                               #markerfacecolor='None',
                               #markeredgecolor=skins['O']['mc'],
                               linestyle=skins['calsource']['ls'],
                               linewidth=skins['calsource']['lw'],
                               #marker=skins['O']['marker'],
                               #markersize=10, 
                               label='Fiber B')                                                   
      ucac2_legend = mlines.Line2D([], [],
                                    color='crimson',
                                    markerfacecolor='None', 
                                    markeredgecolor='crimson', 
                                    linestyle='',
                                    linewidth=0,
                                    marker='o',
                                    markersize=10, label='UCAC2')                         
      
      PM_legend =  mlines.Line2D([], [],color='crimson',
                                 markerfacecolor='crimson',
                                 markeredgecolor='crimson', 
                                 linestyle='-',
                                 linewidth=0.75,
                                 marker='.',
                                 #markersize=10, 
                                 label='PM* (track:$-$%.1f yr)' % (fcm_m.pm_track_time.to(u.yr).value))
                                                          
      ax2.ax.legend(handles=[acq_legend, target_legend, calsource_legend, ucac2_legend, PM_legend],
                 bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                 ncol=6, mode="expand", borderaxespad=0., fontsize=10, borderpad=0.7,
                 handletextpad=0.2, handlelength=2.0)          
      
      
          
# ----------------------------------------------------------------------------------------                    
                    
                    
                    
                    