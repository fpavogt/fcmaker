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
This file contains tools related to the MUSE instrument.
Created October 2017, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''


# ---------------- Define the MUSE field of views -------------------------------------
# From the MUSE manual, March 2014
muse_wfm = np.array([[0.0086683878,   0.0084439565],
                    [359.99128-360., 0.008499508],
                    [359.99172-360.,-0.0082782215],
                    [0.0080579666,  -0.008389310]])
muse_wfm *= 3600 # make it in arcsec

muse_nfm = muse_wfm/8. # For now, just scale the field to 7.5 arcsec ... until I know better


muse_sgs = np.zeros((16,2))
muse_sgs = np.array([[0.012835385,0.0085544861],
                    [0.0085551081,0.012277302],
                    [0.00027758977,0.012],
                    [359.99189-360.,0.012277352],
                    [359.98778-360.,0.0091101444],
                    [359.98739-360.,-0.0091121412],
                    [359.98955-360.,-0.010111818],
                    [359.99417-360.,-0.01050022],
                    [359.99422-360.,-0.011722438],
                    [359.99611-360.,-0.01172232],
                    [359.99622-360.,-0.012166759],
                    [0.004002035,-0.012222326],
                    [0.0040555527,-0.011777885],
                    [0.0060549064,-0.011833571],
                    [0.0060015071,-0.010611345],
                    [0.0090006082,-0.010444971],
                    [0.012444522,-0.0095565611]])
muse_sgs *= 3600 # make it in arcsec

# ----------------------------------------------------------------------------------------

# The default background image from SkyView
bk_image_wfm = 'DSS2 Red'
bk_image_nfm = 'Gaia'

# The radius of the charts
right_radius = 720. # in arcsec

def left_radius(ins_mode):
   ''' A function that sets the defaut size of the left window depending on the mode.'''

   if ins_mode[:3] == 'WFM':
      return 110
   elif ins_mode[:3] == 'NFM':
      return 6
   else:
      raise Exception('Mode undefined.')
      
inner_GS_search = 120. # inner limit to find Guide Stars

# TTS validity area
inner_TTS_search = 1.725/2.*60 # in degree
outer_TTS_search = 3.590/2.*60 # in degree
outer_OATT_search = 3.35 # in arcsec, radius of the on-axis TT star (NFM

# List the supported MUSE observing templates
muse_templates = [# WFM-NOAO
                  'MUSE_wfm-noao_acq_preset',
                  'MUSE_wfm-noao_acq_presetRRM',
                  'MUSE_wfm-noao_acq_movetopixel',
                  'MUSE_wfm-noao_obs_genericoffset',
                  # WFM-AO
                  'MUSE_wfm-ao_acq_movetopixelLGS',
                  'MUSE_wfm-ao_obs_genericoffsetLGS',
                  # WFM_cal
                  'MUSE_wfm_cal_specphot',
                  'MUSE_wfm_cal_astrom',
                  #NFM
                  'MUSE_nfm-ao_acq_LGS',
                  'MUSE_nfm-ao_obs_genericoffsetLGS',
                  ]

# The acquisition parameters that matter for the finding charts. Only one dictionary for
# AO and NOAO - we then update only the parameters that matter.

muse_acq_params = {'ins_mode': None, # Instrument mode
                   'bos_ra': 0,      # blind offset RA
                   'bos_dec': 0,     # blind offset Dec
                   'acq_pa': 0,      # PA
                   'is_gs': False,   # Guide star defined by user ?
                   'gs': None,       # gs coord as SkyCoord
                   'is_tts':False,   # TTS provided ?
                   'tts1': None,     # tts1 coord as SkyCoord
                   'tts2': None,     # tts2 coord as SkyCoord
                   'oatt':None,      # on-axis TTS for NFM
                   'movetopix':True,# whether it is a movetopixel or preset acquisition. Assume movetopix by default, including for "local" files.
                  }
         
muse_sci_params = {'noff': 1,        # Number of offsets
                   'obstype': ['O'], # Observation type (S,O)
                   'coordtype': 'SKY', # Offset type (SKY, DETECTOR)
                   'posang': [0],    # Position angles
                   'off1': [0],      # RA offsets
                   'off2': [0],      # DEC offsets
                   'return': True,   # return to origin ?
                   }
               

# ----------------------------------------------------------------------------------------
def detector_to_sky(dx,dy,pa):
   '''
   Converts MUSE offsets from the DETECTOR reference frame to the SKY reference frame.
   
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
def get_p2fcdata_muse(fc_params, ob, api):
   '''
   Extracts all the important info to build a finding chart from a given MUSE OB from p2.
   
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
      if not(t['templateName'] in muse_templates):
         warnings.warn('Template %s is not supported by fcmaker.' % (t['templateName']))
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
         fc_params['acq'] = copy.deepcopy(muse_acq_params)
         
         tpl,tplVersion = api.getTemplate(obId, t['templateId'])
         
         # If this a preset or movetopix template ?
         if 'preset' in t['templateName']:
            fc_params['acq']['movetopix'] = False
         
         for param in tpl['parameters']:
               
            # Instrument mode
            if param['name'] == 'INS.MODE':
               fc_params['ins_mode'] = param['value'] # The mode is set in the acquisition
            
            # The NFM has got no 'ins_mode' anymore ... deal with it the other way
            if 'nfm' in t['templateName']:
               fc_params['ins_mode'] = 'NFM'
                  
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
            if param['name'] == 'INS.DROT.POSANG':             
               fc_params['acq']['acq_pa'] = param['value']    
         
            # TTSs
            if param['name'] == 'SEQ.AOF.TTS':             
               fc_params['acq']['is_tts'] = param['value']
            if param['name'] == 'TEL.TTS1.ALPHA':             
               tts1_ra = param['value']   
            if param['name'] == 'TEL.TTS1.DELTA':             
               tts1_dec = param['value']
            if param['name'] == 'TEL.TTS2.ALPHA':             
               tts2_ra = param['value']   
            if param['name'] == 'TEL.TTS2.DELTA':             
               tts2_dec = param['value']  
            
            # NFM TTS
            if param['name'] == 'SEQ.NGS.ALPHA':             
               oatt_ra = param['value']   
            if param['name'] == 'SEQ.NGS.DELTA':             
               oatt_dec = param['value']
            
            
         # Store the GS and TTS as SkyCoords
         if fc_params['acq']['is_gs']:
            fc_params['acq']['gs'] = SkyCoord(gs_ra, gs_dec, 
                                       obstime = fcm_m.obsdate, 
                                       #equinox=ob['target']['equinox'], 
                                       frame='icrs',unit=(u.hourangle, u.deg))
                                       
         if fc_params['acq']['is_tts']:
            fc_params['acq']['tts1'] = SkyCoord(tts1_ra,tts1_dec, 
                                         obstime = fcm_m.obsdate, 
                                         #equinox=ob['target']['equinox'], 
                                         frame='icrs',unit=(u.hourangle, u.deg))
            fc_params['acq']['tts2'] = SkyCoord(tts2_ra,tts2_dec, 
                                         obstime = fcm_m.obsdate, 
                                         #equinox=ob['target']['equinox'], 
                                         frame='icrs',unit=(u.hourangle, u.deg))
      
         # Also for the NFM case
         if 'NFM' in fc_params['ins_mode']:
            fc_params['acq']['oatt'] = SkyCoord(oatt_ra,oatt_dec, 
                                         obstime = fcm_m.obsdate, 
                                         #equinox=ob['target']['equinox'], 
                                         frame='icrs',unit=(u.hourangle, u.deg))
      
      
      elif t['type'] in ['science','calib']:
         
         # Ok, I found one more Science template
         fc_params['n_sci'] += 1
         
         temp_name = 'sci%i' %(fc_params['n_sci'])
         
         # Store all the relevant parameters in a dictionary
         fc_params[temp_name] = copy.deepcopy(muse_sci_params)
         
         tpl,tplVersion = api.getTemplate(obId, t['templateId'])
         for param in tpl['parameters']:
               
            if param['name'] == 'SEQ.NOFF':             
               fc_params[temp_name]['noff'] = param['value']
            if param['name'] == 'SEQ.OBSTYPE.LIST':             
               fc_params[temp_name]['obstype'] = param['value']
            if param['name'] == 'SEQ.OFFSET.COORDS':             
               fc_params[temp_name]['coordtype'] = param['value']      
            if param['name'] == 'SEQ.OFFSET.POSANG.LIST':             
               fc_params[temp_name]['posang'] = param['value']  
            if param['name'] == 'SEQ.OFFSET1.LIST':             
               fc_params[temp_name]['off1'] = param['value']
            if param['name'] == 'SEQ.OFFSET2.LIST':             
               fc_params[temp_name]['off2'] = param['value']
            if param['name'] == 'SEQ.RETURN':
               fc_params[temp_name]['return'] = param['value']
   
                          
   return fc_params

# ------------------------------------------------------------------------------
def get_localfcdata_muse(fc_params,inpars):
   '''
   Extracts all the important info to build a finding chart from a given MUSE OB defined
   locally.
   
   Args:
      inpars: A dictionnary containing the OB parameters
   Returns:
       A dictionnary containing the MUSE OB parameters
   ''' 
  
   fc_params['ins_mode'] = inpars['ins_mode']                                     
         
   # Acquisition
   fc_params['acq'] = copy.deepcopy(muse_acq_params)
   fc_params['acq']['is_tts'] = inpars['is_tts']
   # Override the user input if clash with instrument mode
   if 'WFM-NOAO' in inpars['ins_mode']:
      inpars['is_tts'] = False   
   fc_params['acq']['bos_ra'] = inpars['bos_ra']
   fc_params['acq']['bos_dec'] = inpars['bos_dec']
   fc_params['acq']['is_gs'] = inpars['is_gs']
   if inpars['is_gs']:
      fc_params['acq']['gs'] = SkyCoord(inpars['gs_ra'],inpars['gs_dec'], 
                                        frame = 'icrs', 
                                        obstime = Time(fcm_m.obsdate),
                                        equinox = 'J2000')
   if 'WFM-AO' in inpars['ins_mode']:
   
      fc_params['acq']['tts1'] = SkyCoord(inpars['tts1_ra'],inpars['tts1_dec'], 
                                          frame = 'icrs',
                                          obstime = Time(fcm_m.obsdate),
                                          equinox = 'J2000')
      fc_params['acq']['tts2'] = SkyCoord(inpars['tts2_ra'],inpars['tts2_dec'], 
                                          frame='icrs', 
                                          obstime = Time(fcm_m.obsdate),
                                          equinox='J2000')
   
   if 'NFM' in inpars['ins_mode']:
      
      if len(inpars['oatt_ra'])==0 or len(inpars['oatt_dec'])==0:
         raise Exception('Ouch! In NFM, you *must* specify the on-axis tip-tilt (OATT) star coordinates!')
         
      fc_params['acq']['oatt'] = SkyCoord(inpars['oatt_ra'],inpars['oatt_dec'], 
                                          frame='icrs', 
                                          obstime = Time(fcm_m.obsdate),
                                          equinox='J2000')
                                                               
   fc_params['acq']['acq_pa'] = inpars['acq_pa']

   # Observation
   fc_params['n_sci'] = 1
   fc_params['sci1'] = copy.deepcopy(muse_sci_params)
   fc_params['sci1']['noff'] = inpars['noff']
   fc_params['sci1']['obstype'] = [i for i in inpars['obstype'][0].split(' ')]
   fc_params['sci1']['posang'] = [float(i) for i in str(inpars['posang'][0]).split(' ')]
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
   counter = 1
   acq_field =  fcm_t.offset_coord(target, 
                                   delta_ra=fc_params['acq']['bos_ra']*u.arcsec, 
                                   delta_dec=fc_params['acq']['bos_dec']*u.arcsec)
   
   # Start filling the dictionary of fields
   fields[1] = [fc_params['inst'], # Instrument first
                fc_params['ins_mode'], # Instrument mode
                acq_field, # Field central coordinates
                fc_params['acq']['acq_pa'], # Position Angle
                'Acq', # Nature of the field ('Acq', 'Target', 'O', 'S')
                ]
      
   # Store the PA of the acquisition
   delta_pa = fc_params['acq']['acq_pa']
      
   # If we are in MUSE WFM AO, the system will close the loop on the target. 
   # So here, I add it as one field, EVEN IF the user then applies offsets in the 
   # observation templates. This allows to check that the TTS will be ok, when we will 
   # first close the loop.
   counter+=1
   fields[counter] = [fc_params['inst'], # Instrument first
                      fc_params['ins_mode'], # Instrument mode
                      target, # Field central coordinates
                      fc_params['acq']['acq_pa'], # Position Angle
                      'Target', # Nature of the field ('Acq', 'Target', 'O', 'S')
                      ]  

   # Now, loop through all the Science templates I have
   for n in range(fc_params['n_sci']):

      # What is the name of this specific template ?
      temp_name = 'sci%i' %(n+1)
      
      # Ok, what are the offset sequences
      off1 = fc_params[temp_name]['off1']
      off2 = fc_params[temp_name]['off2']
      posang = fc_params[temp_name]['posang']
      obstype = fc_params[temp_name]['obstype']
      coordtype = fc_params[temp_name]['coordtype']
      
      # If any sequence is too short, loop it as required
      all_offs = [off1, off2, posang, obstype]
      for (s,seq) in enumerate(all_offs):

         if len(seq) < fc_params[temp_name]['noff']:
            seq = seq * (fc_params[temp_name]['noff']//len(seq)+1)
            all_offs[s] = seq[:fc_params[temp_name]['noff']]
  
         #elif len(seq) > fc_params[temp_name]['noff']:
         #   warnings.warn('Offset sequence larger than "NOFF". Ignoring extra bits ...')
      
      [off1, off2, posang, obstype] = all_offs

      # Then build the OB sequence
      for o in range(fc_params[temp_name]['noff']): 
         
         counter +=1
         
         
         if coordtype == 'SKY':
            # Sum the RA, Dec offsets
            delta_ra += off1[o]
            delta_dec += off2[o]
            # Sum the PA offsets
            delta_pa += posang[o]
            
         elif coordtype == 'DETECTOR':
            # Transform the dx, dy values in dra, ddec, given the current PA
            (this_dra,this_ddec) = detector_to_sky(off1[o],off2[o],delta_pa)
            
            # Then sum these and apply any PA shift
            delta_ra += this_dra
            delta_dec += this_ddec
            # Sum the PA offsets
            delta_pa += posang[o]
         
         else:
            raise Exception('Ouch! coordtype unknwown: %s' % (coordtype))
           
         # Create the field entry
         fields[counter] = [fc_params['inst'],
                            fc_params['ins_mode'],
                            fcm_t.offset_coord(target, 
                                               delta_ra = delta_ra*u.arcsec,
                                               delta_dec = delta_dec*u.arcsec,),
                            delta_pa % 360,
                            obstype[o]
                           ]
   
      # Ok, I return to the origin after the end of the OB, so let's reset all the offsets
      # to match the acquisition.
      if fc_params[temp_name]['return']:
         delta_ra = 0
         delta_dec = 0
         delta_pa = fc_params['acq']['acq_pa']
   
   return fields
   
# ------------------------------------------------------------------------------
def get_polygon(central_coord, pa, mode):
   '''
   Given the central location and position of a field, build a polygon to feed to 
   matplotlib down the line.
   
   Args:
      central_coord: an astropy.SkyCoord entry with the center of the MUSE field
      pa: the position angle of the field (p2 convention)
      mode: the instrument mode, either 'WFM' or 'NFM'
   Returns:
      A list of 2-D coordinates for each corner of the field-of-view.
   '''
       
   # Very well, having done this, I will also build a list of polygons for each field, 
   # that I can feed directly to matplotlib
   polygon = [] #it's a list
   
   
   if mode == 'WFM':
      corners = copy.deepcopy(muse_wfm)
   elif mode == 'NFM':
      corners = copy.deepcopy(muse_nfm)
   else:
      raise Exception('Mode unknown.')
       
   # Rotate the different corners by the correct amount
   for (j,corner) in enumerate(corners):    
      corners[j] = np.dot(corners[j],fcm_t.myrotmatrix(pa+180))
        
   # Converts from arc sec to degrees/hours
   corners = corners /3600.
    
   # Finds the field center
   thisalpha = central_coord.ra
   thisdelta = central_coord.dec
    
   # Accounts for the declination in the RA size of the field
   corners[:,0] = corners[:,0]/np.cos(thisdelta.radian)

   # Store the field corners
   for (j,corner) in enumerate(corners):
      polygon.append(corners[j] + [thisalpha.deg,thisdelta.deg])
    
   polygon = [np.array(polygon)]
                     
   return polygon
   
# ----------------------------------------------------------------------------------------
def plot_field(ax1, ax2, fc_params, field):
   '''
   The specific MUSE function that draws a specific observation field.
   
   Args:
      ax1,ax2: the left and right plots 'axis'
      fc_params: the dictionnary of parameters for the OB
      field: the specific field list of parameters
   '''
   
   skins = {'Acq':{'c':'darkorchid', 'lwm':2, 'zorder':10, 'marker':fcm_t.crosshair(pa=45), 
                   'lw':1.5, 'ms':500, 'ls':'-', 'mc':'darkorchid'},
            'Target': {'c':'None', 'lwm':1, 'zorder':5, 'marker':'D', 
                       'lw':1.5, 'ms':250, 'ls':'None', 'mc':'darkorange'},
            'O':{'c':'royalblue', 'lwm':1, 'zorder':5, 'marker':'o', 
                 'lw':1.5, 'ms':50, 'ls':':', 'mc':'royalblue'},
            'S':{'c':'darkcyan', 'lwm':1, 'zorder':5, 'marker':'s', 
                 'lw':1.5, 'ms':50, 'ls':':', 'mc':'darkcyan'},
            }
   
   
   this_coords = [field[2].ra, field[2].dec]
   
   for ax in [ax1,ax2]:
      # Show markers for NFM always, or only for ax1
      if (ax in [ax1]) or (fc_params['ins_mode'] == 'NFM'):
         # Do not show the marker for "preset" acquisitions
         if (not(field[4]) is 'Acq') or fc_params['acq']['movetopix']:
      
            # Show the center of the field, without obstructing it to see the blind offset star
            ax.show_markers(this_coords[0].deg, this_coords[1].deg, 
                            marker=skins[field[4]]['marker'],
                            edgecolor=skins[field[4]]['mc'],
                            s=skins[field[4]]['ms'], 
                            linewidth=skins[field[4]]['lwm'], 
                            zorder=skins[field[4]]['zorder'],
                            ) 
   
      # instrument footprint, except for the "target", where the AO loops are closed, but
      # no data is taken
      #if field[4] != 'Target':
      if (ax in [ax1]) or (fc_params['ins_mode'] != 'NFM'):
         ax.show_polygons(get_polygon(field[2],field[3],field[1][:3]), 
                           edgecolor = skins[field[4]]['c'],
                           linewidth = skins[field[4]]['lw'],
                           zorder = skins[field[4]]['zorder'],
                           linestyle = skins[field[4]]['ls'],)
            
      # the GS validity area
      ax.show_circles([this_coords[0].deg],
                      [this_coords[1].deg],
                      [((fcm_m.outer_GS_Nas*u.arcsec).to(u.degree)).value,],
                       color='k', lw=0.5, 
                       zorder=skins[field[4]]['zorder']
                     )           
                                    
      # TTS, GS footprint for the Science only - no TT star used for the sky fields, 
      # Only show the TT area for AO mode
      if (field[4] in ['Target','O']) and ('WFM-AO' in fc_params['ins_mode']):
            
         ax.show_circles([this_coords[0].deg,this_coords[0].deg],
                         [this_coords[1].deg,this_coords[1].deg],
                         [((inner_TTS_search*u.arcsec).to(u.degree)).value,
                          ((outer_TTS_search*u.arcsec).to(u.degree)).value,
                         ],
                         color='k', lw=0.5, linestyle= '-')  
                         
         # TODO: here, also show the valid area of the TTS      
                  
      # Show the NFM TT star
      if 'NFM' in fc_params['ins_mode']:
         tts = fc_params['acq']['oatt']  
            
         # Skip if it is not defined (0,0 is the default)
         if tts.ra == 0 and tts.dec == 0:
            continue
                          
         ax.show_markers(tts.ra, tts.dec, marker=fcm_t.crosshair(pa=0), edgecolor='tomato',
                           s=500, linewidth=1.5)
         
         # Only add their name to the zoom-in plot
         if fcm_m.fcm_usetex:
            lab = r'\textbf{TT}'
         else:
            lab = r'TT'
         if ax == ax1 and (field[4] == 'Acq'):
            ax.add_label(tts.ra.deg,tts.dec.deg+0.75/3600,lab, 
                         verticalalignment='center', 
                         horizontalalignment='center',size=10,color='k',
                         bbox=dict(facecolor='w',ec='k', alpha=0.6)) 
            
         # Check if the TTS is compatible with this field. If not, flag it as such !
         # Make sure this is a 'O' or 'Target' field as well ...
         if ax ==ax1 and \
            (field[4] in ['Target','O']) and \
            (tts.separation(field[2]) > (outer_OATT_search*u.arcsec)):
                
            if fcm_m.fcm_usetex:
               lab = r'\textbf{!}'
            else:
               lab = '!'
                  
            ax.add_label(tts.ra.deg,tts.dec.deg,lab, 
                            verticalalignment='center', 
                            horizontalalignment='center',size=12,color='tomato',
                            bbox=dict(boxstyle="circle,pad=0.17", facecolor='w',ec='tomato', alpha=1))
         
         # Show the OATT validity area                
         if (field[4] in ['Target','O']):
            
            ax.show_circles([this_coords[0].deg],
                            [this_coords[1].deg],
                            [((outer_OATT_search*u.arcsec).to(u.degree)).value,],
                            color='k', lw=0.5, linestyle= '-')                 
                           
                           
      # Show the WFM-AO TT stars 
      if fc_params['acq']['is_tts'] and ('WFM-AO' in fc_params['ins_mode']):
         for (t,tts) in enumerate([fc_params['acq']['tts1'],fc_params['acq']['tts2']]):  
            
            # Skip if it is not defined (0,0 is the default)
            if tts.ra == 0 and tts.dec == 0:
               continue
                          
            ax.show_markers(tts.ra, tts.dec, marker=fcm_t.crosshair(pa=0), edgecolor='tomato',
                               s=500, linewidth=1.5)
   
            # Only add their name to the zoom-in plot
            if fcm_m.fcm_usetex:
               lab = r'\textbf{TTS%i}'%(t+1)
            else:
               lab = r'TTS%i'%(t+1)
            if ax == ax1 and (field[4] == 'Acq'):
               ax.add_label(tts.ra.deg,tts.dec.deg+15./3600,lab, 
                            verticalalignment='center', 
                            horizontalalignment='center',size=12,color='k',
                            bbox=dict(facecolor='w',ec='k', alpha=0.6)) 
            
            # Check if the TTS is compatible with this field. If not, flag it as such !
            # Make sure this is a 'O' or 'Target' field as well ...
            if ax ==ax1 and \
               (field[4] in ['Target','O']) and \
               ((tts.separation(field[2]) > (outer_TTS_search*u.arcsec)) or \
                (tts.separation(field[2]) < (inner_TTS_search*u.arcsec))):
                
               if fcm_m.fcm_usetex:
                  lab = r'\textbf{!}'
               else:
                  lab = '!'
                  
               ax.add_label(tts.ra.deg,tts.dec.deg,lab, 
                            verticalalignment='center', 
                            horizontalalignment='center',size=12,color='tomato',
                            bbox=dict(boxstyle="circle,pad=0.17", facecolor='w',ec='tomato', alpha=1)) 
      
          
                   
      # Show the Guide Star (if there is an acq template present)
      if (fc_params['acq']['is_gs']) and (field[4] == 'Acq'):
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
   O_legend = mlines.Line2D([], [], 
                               color=skins['O']['c'],
                               markerfacecolor='None',
                               markeredgecolor=skins['O']['mc'],
                               linestyle=skins['O']['ls'],
                               linewidth=skins['O']['lw'],
                               marker=skins['O']['marker'],
                               markersize=10, label='O')                                                           
   S_legend = mlines.Line2D([], [], 
                               color=skins['S']['c'],
                               markerfacecolor='None',
                               markeredgecolor=skins['S']['mc'], 
                               linestyle=skins['S']['ls'],
                               linewidth=skins['S']['lw'],
                               marker=skins['S']['marker'],
                               markersize=10, label='S') 
                               
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
                                                          
   ax2._ax1.legend(handles=[acq_legend, target_legend,O_legend,S_legend,ucac2_legend, PM_legend],
                 bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                 ncol=6, mode="expand", borderaxespad=0., fontsize=10, borderpad=0.7,
                 handletextpad=0.2, handlelength=2.0)          
                    
# ----------------------------------------------------------------------------------------                    
                    
                    
                    
                    