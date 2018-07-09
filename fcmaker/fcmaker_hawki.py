# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------------------

import numpy as np
from astropy.coordinates.sky_coordinate import SkyCoord
import astropy.units as u
from astropy.time import Time
import warnings
import copy


import matplotlib.lines as mlines

from . import fcmaker_tools as fcm_t
from . import fcmaker_metadata as fcm_m

'''
fcmaker: a Python module to automatically create finding charts for ESO OBs in p2.\n
Copyright (C) 2017,  F.P.A. Vogt
--- oOo ---
This file contains tools related to the HAWKI instrument.
Created October 2017, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''


# ---------------- Define the HAWKI field of view -------------------------------------
# From the HAWKI manual, June 2017
hawki_pix = 0.106 # size of a HAWKI pixel
                     
# The rotation of each quadrant [Q1,Q2,Q3,Q4], in degrees                     
Q_rot = np.array([0, 0.13, 0.03, 0.04]) 

# Spatial offset of quadrants from telescope pointing and rotation center [Q1, Q2, Q3, Q4]
Q_off =[np.array([-2163., -2164.]),
        np.array([37.5, -2165.5]),
        np.array([42.,28.]),
        np.array([-2158.,25.5])] # in pixels                              

# ----------------------------------------------------------------------------------------

# The default background image from SkyView
bk_image = '2MASS-H'

# The radius of the charts
right_radius = 720. # in arcsec
left_radius = 310. # in arcsec
inner_GS_search = 4. * 60. # inner limit to find Guide Stars

#: List all the supported HAWKI observing templates (incl. Fast Phot)
hawki_templates = [# NOAO
                  'HAWKI_img_acq_Preset',
                  'HAWKI_img_acq_PresetRRM',
                  'HAWKI_img_acq_MoveToPixel',
                  'HAWKI_img_acq_FastPhot',
                  'HAWKI_img_acq_FastPhotRRM',
                  'HAWKI_img_obs_GenericOffset',
                  'HAWKI_img_obs_AutoJitter',
                  # AO
                  'HAWKI_img_acq_LGS_Preset',
                  'HAWKI_img_acq_LGS_PresetRRM',
                  'HAWKI_img_acq_LGS_MoveToPixel',
                  'HAWKI_img_acq_LGS_FastPhot',
                  # cal
                  ]
                  
# List of FastPhot templates only                   
# Not needed as long as the string 'FastPhot' is common to all these.
#hawki_fph_templates = ['HAWKI_img_acq_FastPhot',
#                       'HAWKI_img_acq_FastPhotRRM',]

# The acquisition parameters that matter for the finding charts. Only one dictionary for
# AO and NOAO - we then update only the parameters that matter.

hawki_acq_params = {'filter': None,   # Filter name
                    'bos_ra': 0,      # Alpha offset for the target
                    'bos_dec': 0,     # Delta offset for the target
                    'acq_pa': 0,      # PA
                    'is_gs': False,   # Guide star defined by user ?
                    'gs': None,       # gs coord as SkyCoord
                    'is_fph': False,  # A tag to keep track of FastPhot OBs
                    'nx':128,        
                    'ny':2048,
                    'startx':1,
                    'starty':1,
                  }
         
hawki_sci_params = {'noff': 0,        # Number of offsets
                    'filter': None, # Filter name
                    'obstype': ['O'], # Observation type (S,O)
                    'coordtype': 'SKY', # Offset type (SKY, DETECTOR)
                    'off1': [0],      # RA offsets
                    'off2': [0],      # DEC offsets
                    'return': True,   # return to origin ?
                    'jitter_box': 0,
                    'ra_sky': None,   # Fixed offset to sky
                    'dec_sky': None,  # Fixed offset to sky
                    # The fast phot aspects below can never be set in the SCIENCE templates
                    # but I include them for consistency in sub-functions down the line.
                    'is_fph': False,  # A tag to keep track of FastPhot OBs
                    'nx':128,        
                    'ny':2048,
                    'startx':1,
                    'starty':1,
                   }
               

# ----------------------------------------------------------------------------------------
def detector_to_sky(dx,dy,pa):
   '''
   Converts HAWKI offsets from the DETECTOR reference frame to the SKY reference frame.
   
   Args:
      dx (float): offset in arcsec
      dy (float): offset in arcsec
      pa (float): current position angle in degrees measured East-of-North (MUSE convention)
        
   Returns:
      tuple of floats: the (dra,ddec) offsets in SKY convention.
      
   Note:
      According to the MUSE User manual, "when in DETECTOR framework, if both telescope
      offsets and PA offsets are provided, the telescope offsets are sent (first) with 
      respect to the current PA of the instrument, before applying the PA offset."
   '''
   
   # Rotate the offsets in the Ra-Dec frame
   ddec = np.sin(np.radians(pa)) * dx + np.cos(np.radians(pa)) * dy
   dra = np.cos(np.radians(pa)) * dx - np.sin(np.radians(pa)) * dy
   
   # Flip the RA axis to increase along East
   dra *= -1.
   
   return (dra, ddec)

# ----------------------------------------------------------------------------------------
def get_p2fcdata_hawki(fc_params, ob, api):
   '''
   Extracts all the important info to build a finding chart from a given HAWKI OB from p2.
   
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
      if not(t['templateName'] in hawki_templates):
         warnings.warn('Template %s is not supported by fcmaker. Ignoring it.' % (t['templateName']))
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
         fc_params['acq'] = copy.deepcopy(hawki_acq_params)
         
         # If this is a FastPhot OB, take note of it.
         #if t['templateName'] in hawki_fph_templates:
         if 'FastPhot' in t['templateName']:
            fc_params['acq']['is_fph'] = True
         
         tpl,tplVersion = api.getTemplate(obId, t['templateId'])
         for param in tpl['parameters']:
               
            # Instrument mode
            if param['name'] == 'INS.FILT.NAME':
               fc_params['acq']['filter'] = param['value'] # The filter set in the acquisition
                  
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
               
            # Position Angle
            if param['name'] == 'TEL.ROT.OFFANGLE':             
               fc_params['acq']['acq_pa'] = param['value']    
            
            # The FastPhot parameters
            if param['name'] == 'SEQ.WIN.NX':             
               fc_params['acq']['nx'] = param['value'] 
               
               # For now, take the easy way out if nx != 128
               if fc_params['acq']['nx'] != 128:
                  raise Exception('Ouch! In Fast Phot, nx=128 is required!') 
                  
            if param['name'] == 'SEQ.WIN.NY':             
               fc_params['acq']['ny'] = param['value'] 
            if param['name'] == 'SEQ.WIN.STARTX':             
               fc_params['acq']['startx'] = param['value'] 
            if param['name'] == 'SEQ.WIN.STARTY':             
               fc_params['acq']['starty'] = param['value'] 
              
         
         # Store the GS and TTS as SkyCoords
         if fc_params['acq']['is_gs']:
            fc_params['acq']['gs'] = SkyCoord(gs_ra, gs_dec, 
                                       obstime = fcm_m.obsdate, 
                                       #equinox=ob['target']['equinox'], 
                                       frame='icrs',unit=(u.hourangle, u.deg))
                                       
      
      elif t['type'] in ['science']:
         
         # Ok, I found one more Science template
         fc_params['n_sci'] += 1
         
         temp_name = 'sci%i' %(fc_params['n_sci'])
         
         # Store all the relevant parameters in a dictionary
         fc_params[temp_name] = copy.deepcopy(hawki_sci_params)
         
         tpl,tplVersion = api.getTemplate(obId, t['templateId'])
         for param in tpl['parameters']:
               
            if param['name'] == 'SEQ.NOFFSET':             
               fc_params[temp_name]['noff'] = param['value']
            if param['name'] == 'SEQ.OBSTYPE.LIST':             
               fc_params[temp_name]['obstype'] = param['value']
            if param['name'] == 'SEQ.OFFSET.COORDS':             
               fc_params[temp_name]['coordtype'] = param['value']       
            if param['name'] == 'SEQ.OFFSET1.LIST':             
               fc_params[temp_name]['off1'] = param['value']
            if param['name'] == 'SEQ.OFFSET2.LIST':             
               fc_params[temp_name]['off2'] = param['value']
            if param['name'] == 'SEQ.RETURN':
               fc_params[temp_name]['return'] = param['value']
            if param['name'] == 'INS.FILT.NAME':
               fc_params[temp_name]['filter'] = param['value']
   
                          
   return fc_params

# ------------------------------------------------------------------------------
def get_localfcdata_hawki(fc_params, inpars):
   '''
   Extracts all the important info to build a finding chart from a given MUSE OB defined
   locally.
   
   Args:
      inpars: A dictionnary containing the OB parameters
   Returns:
       A dictionnary containing the MUSE OB parameters
   ''' 
   
   # Acquisition
   fc_params['acq'] = copy.deepcopy(hawki_acq_params)
   fc_params['acq']['filter'] = inpars['acq_filter'] 
   fc_params['acq']['bos_ra'] = inpars['bos_ra']
   fc_params['acq']['bos_dec'] = inpars['bos_dec']
   fc_params['acq']['is_gs'] = inpars['is_gs']
   fc_params['acq']['gs'] = SkyCoord(inpars['gs_ra'],inpars['gs_dec'], 
                                     frame = 'icrs', 
                                     obstime = Time(fcm_m.obsdate),
                                     equinox = 'J2000')
   fc_params['acq']['acq_pa'] = inpars['acq_pa']
   fc_params['acq']['is_fph'] = inpars['is_fph']
   
   if inpars['is_fph']:
      fc_params['acq']['nx'] = inpars['nx']
   
      # For now, take the easy way out if nx != 128
      if fc_params['acq']['nx'] != 128:
         raise Exception('Ouch! In Fast Phot, nx=128 is required!') 
      
      fc_params['acq']['ny'] = inpars['ny']
      fc_params['acq']['startx'] = inpars['startx']
      fc_params['acq']['starty'] = inpars['starty']
   
   # Observation
   fc_params['n_sci'] = 1
   fc_params['sci1'] = copy.deepcopy(hawki_sci_params)
   fc_params['sci1']['noff'] = inpars['noff']
   fc_params['sci1']['obstype'] = [i for i in inpars['obstype'][0].split(' ')]
   fc_params['sci1']['off1'] = [float(i) for i in str(inpars['off1'][0]).split(' ')]
   fc_params['sci1']['off2'] = [float(i) for i in str(inpars['off2'][0]).split(' ')]
   fc_params['sci1']['return'] = inpars['return']
   fc_params['sci1']['coordtype'] = inpars['coordtype']

   return fc_params

# ------------------------------------------------------------------------------
def get_fields_dict(fc_params):
   '''
   Create a dictionnary with the basic info required to draw the fields, given 
   certain OB parameters.
   
   Args:
      fc_params: A dictionnary containing the MUSE OB parameters
   Returns:
       A dictionary containing the plotting parameters for each field.
   
   '''
   
   # Some values to keep track of the ra, dec and pa offsets over the course of the OB.
   delta_ra = 0
   delta_dec = 0
   delta_pa = 0

   # Define the center of the finding chart to be on the acquisition field.
   target = fc_params['target']  
   
   # Observing fields : build a dictionary for all of them
   fields = {}

   # First, the acquisiton frame 
   # Am I counting the offsets in the correct way ?
   counter = 1
   
   delta_ra += fc_params['acq']['bos_ra']
   delta_dec += fc_params['acq']['bos_dec']
   
   acq_field =  fcm_t.offset_coord(target, 
                                   delta_ra = delta_ra * u.arcsec, 
                                   delta_dec = delta_dec * u.arcsec)
   
   # Start filling the dictionary of fields
   fields[1] = [fc_params['inst'], # Instrument first
                fc_params['acq']['filter'], # Filter
                acq_field, # Field central coordinates
                fc_params['acq']['acq_pa'], # Position Angle
                'Acq', # Nature of the field ('Acq', 'Target', 'O', 'S')
                fc_params['acq']['is_fph'],
                fc_params['acq']['nx'],
                fc_params['acq']['ny'],
                fc_params['acq']['startx'],
                fc_params['acq']['starty'],
                ]
      
   # Store the PA of the acquisition
   delta_pa = fc_params['acq']['acq_pa']
   
   '''  
   # For information, also show the target field.
   counter+=1
   fields[counter] = [fc_params['inst'], # Instrument first
                      fc_params['acq']['filter'], # Filter
                      target, # Field central coordinates
                      fc_params['acq']['acq_pa'], # Position Angle
                      'Target', # Nature of the field ('Acq', 'Target', 'O', 'S')
                      ]  
   '''
   
   # Now, loop through all the Science templates I have
   for n in range(fc_params['n_sci']):

      # What is the name of this specific template ?
      temp_name = 'sci%i' %(n+1)
      
      # Ok, what are the offset sequences
      off1 = fc_params[temp_name]['off1']
      off2 = fc_params[temp_name]['off2']
      obstype = fc_params[temp_name]['obstype']
      coordtype = fc_params[temp_name]['coordtype']
      
      # If any sequence is too short, loop it as required
      all_offs = [off1, off2, obstype]
      for (s,seq) in enumerate(all_offs):

         if len(seq) < fc_params[temp_name]['noff']:
            seq = seq * (fc_params[temp_name]['noff']//len(seq)+1)
            all_offs[s] = seq[:fc_params[temp_name]['noff']]
  
         #elif len(seq) > fc_params[temp_name]['noff']:
         #   warnings.warn('Offset sequence larger than "NOFF". Ignoring extra bits ...')
      
      [off1, off2, obstype] = all_offs

      # Then build the OB sequence
      for o in range(fc_params[temp_name]['noff']): 
         
         counter +=1
         
         if coordtype == 'SKY':
            # Sum the RA, Dec offsets
            delta_ra += off1[o]
            delta_dec += off2[o]
            # Sum the PA offsets
            #delta_pa += posang[o]
            
         elif coordtype == 'DETECTOR':
            # Transform the dx, dy values in dra, ddec, given the current PA
            (this_dra,this_ddec) = detector_to_sky(off1[o],off2[o],delta_pa)
            
            # Then sum these and apply any PA shift
            delta_ra += this_dra
            delta_dec += this_ddec
            # Sum the PA offsets
            #delta_pa += posang[o]
         
         else:
            raise Exception('Ouch! coordtype unknwown: %s' % (coordtype))
            
         # Create the field entry - am I getting the sign of the offsets right ?
         fields[counter] = [fc_params['inst'],
                            fc_params[temp_name]['filter'],
                            fcm_t.offset_coord(target, 
                                               delta_ra = delta_ra*u.arcsec,
                                               delta_dec = delta_dec*u.arcsec,),
                            delta_pa % 360,
                            obstype[o],
                            fc_params[temp_name]['is_fph'],
                            fc_params[temp_name]['nx'],
                            fc_params[temp_name]['ny'],
                            fc_params[temp_name]['startx'],
                            fc_params[temp_name]['starty'],
                           ]
   
      # Ok, I return to the origin after the end of the OB, so let's reset all the offsets
      # to match the acquisition.
      if fc_params[temp_name]['return']:
         delta_ra = 0
         delta_dec = 0
         delta_pa = fc_params['acq']['acq_pa']
   
   return fields
   
# ------------------------------------------------------------------------------
def get_polygon(central_coord, pa, nx, ny, startx, starty):
   '''
   Given the central location and position of a field, build a polygon to feed to 
   matplotlib down the line.
   
   Args:
      central_coord: an astropy.SkyCoord entry with the center of the MUSE field
      pa: the position angle of the field (p2 convention)
      nx,ny, startx, starty: the Fast Phot windowing parameters
   Returns:
      A list of 2-D coordinates for each corner of the field-of-view.
   '''
       
   # Very well, having done this, I will also build a list of polygons for each field, 
   # that I can feed directly to matplotlib
   polygon = [] #it's a list

   # Loop through all the quadrants
   for quad in range(4):
      
      sub_polygon = []
      
      # Deal with possible Fast Phot window, which is mirrored betwwen quad 1,2 and 3,4
      if quad in [0,1]:
         corners = np.array([[startx-1., starty-1.], 
                             [startx-1., starty-1. + ny-1.], 
                             [startx-1.+16*nx-1., starty-1. + ny-1.],
                             [startx-1.+16*nx-1., starty-1.]])
      elif quad in [2,3]:
         corners = np.array([[startx-1., 2047 - (starty-1.)], 
                             [startx-1., 2047 - (starty-1. + ny-1.)], 
                             [startx-1.+16*nx-1., 2047 - (starty-1. + ny-1.)],
                             [startx-1.+16*nx-1., 2047 - (starty-1.)]])
    
      # Rotate the different corners by the correct amount
      for (j,corner) in enumerate(corners):    
         corners[j] = np.dot(corners[j],fcm_t.myrotmatrix(Q_rot[quad]))
        
      # Shift the array to the correct position with respect to the telescope pointing center
      corners += Q_off[quad]
      
      # Now rotate this again, this time because of the PA
      for (j,corner) in enumerate(corners):    
         # Update v0.2.00 FPAV: fix HAWKI rotation angle
         corners[j] = np.dot(corners[j],fcm_t.myrotmatrix(pa))
      
      # Converts from arc sec to degrees
      corners = corners * hawki_pix/3600.
    
      # Finds the field center
      thisalpha = central_coord.ra
      thisdelta = central_coord.dec
    
      # Accounts for the declination in the RA size of the field
      corners[:,0] = corners[:,0]/np.cos(thisdelta.radian)

      # Store the field corners
      for (j,corner) in enumerate(corners):
         sub_polygon.append(corners[j] + [thisalpha.deg,thisdelta.deg])
    
      polygon = polygon + [np.array(sub_polygon)]
                     
   return polygon
   
# ----------------------------------------------------------------------------------------
def plot_field(ax1, ax2, fc_params, field):
   '''
   The specific HAWKI function that draws a specific observation field.
   
   Args:
      ax1,ax2: the left and right plots 'axis'
      fc_params: the dictionnary of parameters for the OB
      field: the specific field list of parameters
   '''
   
   skins = {'Acq':{'c':'darkorchid', 'lw':2, 'lw_fph':1, 'zorder':10, 'ls':'-'},
            'Target': {'c':'None', 'lwm':1, 'zorder':5, 'marker':'D', 
                       'lw':1.5, 'ms':250, 'ls':'None', 'mc':'darkorange'},
            'O':{'c':'royalblue', 'lw':1, 'zorder':5, 'ls':'--'},
            'S':{'c':'darkcyan', 'lw':1, 'zorder':5, 'ls':'--'},
            }
   
   
   this_coords = [field[2].ra, field[2].dec]
   
   '''
   # Show the center of the field, without obstructing it to see the blind offset star
   ax1.show_markers(this_coords[0].deg, this_coords[1].deg, 
                    marker=skins[field[4]]['marker'],
                    edgecolor=skins[field[4]]['mc'],
                    s=skins[field[4]]['ms'], 
                    linewidth=skins[field[4]]['lwm'], 
                    zorder=skins[field[4]]['zorder'],
                    ) 
   '''
   
   for ax in [ax1,ax2]:
      # full instrument footprint
      ax.show_polygons(get_polygon(field[2],field[3],128.,2048.,1.,1.), 
                       edgecolor = skins[field[4]]['c'],
                       linewidth = skins[field[4]]['lw'],
                       zorder = skins[field[4]]['zorder'],
                       linestyle = skins[field[4]]['ls'],)
      
      # For the acquisition only, show the Quadrant name
      if (field[4] == 'Acq') and (ax == ax2):
      
         xys = np.array([[120.,-120.], [-120.,-120.],
                         [-120.,120.], [120.,120.],])
           
         for (j,xy) in enumerate(xys):    
            xys[j] = np.dot(xys[j],fcm_t.myrotmatrix(fc_params['acq']['acq_pa']))
            
            label_coord = fcm_t.offset_coord(field[2],
                                             delta_ra = xys[j][0]*u.arcsec,
                                             delta_dec = xys[j][1]*u.arcsec)
            
            if fcm_m.fcm_usetex:
               lab = r'\textbf{Q%i}' % (j+1)
            else: 
               lab = r'Q%i' % (j+1)
            
            ax.add_label(label_coord.ra.deg,label_coord.dec.deg,lab, 
                         verticalalignment='center', 
                         horizontalalignment='center',size=12,color=skins[field[4]]['c'],
                         bbox=dict(facecolor='w',ec='none', alpha=0.6), zorder=2
                         ) 
         
      
      # IS this a Fast Phot acquisition ?
      if field[5]:
         # Very well, let's draw the Fast Phot windows
         ax.show_polygons(get_polygon(field[2],field[3],field[6],field[7],field[8],field[9]), 
                          edgecolor = skins[field[4]]['c'],
                          linewidth = skins[field[4]]['lw_fph'],
                          zorder = skins[field[4]]['zorder'],
                          linestyle = skins[field[4]]['ls'],)
      
      
           
      # the GS validity area
      ax.show_circles([this_coords[0].deg],
                      [this_coords[1].deg],
                      [((fcm_m.outer_GS_Nas*u.arcsec).to(u.degree)).value,],
                       color='k', lw=0.5, 
                       zorder=skins[field[4]]['zorder']
                     ) 
                    
      # Show the Guide Star (if there is an acq template present)
      if fc_params['acq']['is_gs']:
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
                                  #markerfacecolor='None',
                                  #markeredgecolor=skins['Acq']['mc'], 
                                  #marker=skins['Acq']['marker'],
                                  linestyle=skins['Acq']['ls'],
                                  linewidth=skins['Acq']['lw'],
                                  #markersize=10, 
                                  label='Acq.')
      '''                            
      target_legend = mlines.Line2D([], [],
                                    color=skins['Target']['c'],
                                    markerfacecolor='None', 
                                    markeredgecolor=skins['Target']['mc'], 
                                    linestyle=skins['Target']['ls'],
                                    linewidth=skins['Target']['lw'],
                                    marker=skins['Target']['marker'],
                                    markersize=10, label='Target')
      '''
      O_legend = mlines.Line2D([], [], 
                               color=skins['O']['c'],
                               #markerfacecolor='None',
                               #markeredgecolor=skins['O']['mc'],
                               linestyle=skins['O']['ls'],
                               linewidth=skins['O']['lw'],
                               #marker=skins['O']['marker'],
                               #markersize=10, 
                               label='O') 
                                                                                         
      S_legend = mlines.Line2D([], [], 
                               color=skins['S']['c'],
                               #markerfacecolor='None',
                               #markeredgecolor=skins['S']['mc'], 
                               linestyle=skins['S']['ls'],
                               linewidth=skins['S']['lw'],
                               #marker=skins['S']['marker'],
                               #markersize=10, 
                               label='S') 
      
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
                                                           
      ax2._ax1.legend(handles=[acq_legend,O_legend,S_legend, ucac2_legend, PM_legend],
                 bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                 ncol=5, mode="expand", borderaxespad=0., fontsize=10, borderpad=0.7,
                 handletextpad=0.2, handlelength=2.0)          
                    
# ----------------------------------------------------------------------------------------                    
                    
                    
                    
                    