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
Copyright (C) 2017,  F.P.A. Vogt
--- oOo ---
This file contains tools related to the MUSE instrument.
Created October 2017, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''


# ---------------- Define the XSHOOTER field of views -------------------------------------
# From the XSHOOTER manual

xshooter_slt_length = 11*u.arcsec

xshooter_ifu = [1.8*u.arcsec,4*u.arcsec]

xshooter_acqcam = [1.5*u.arcmin, 1.5*u.arcmin]

# ----------------------------------------------------------------------------------------

# The default background image from SkyView
bk_image = 'DSS2 Red'

# The radius of the charts
right_radius = 540. # in arcsec

left_radius = 65. # in arcsec
      
inner_GS_search = 120. # inner limit to find Guide Stars # TO BE VALIDATED!
# List the supported MUSE observing templates
xshooter_templates = [# SLT
                      'XSHOOTER_slt_acq',
                      'XSHOOTER_slt_acq_rrm',
                      'XSHOOTER_slt_obs_Stare',
                      'XSHOOTER_slt_obs_StareSynchro',
                      'XSHOOTER_slt_obs_AutoNodOnSlit',
                      'XSHOOTER_slt_obs_FixedSkyOffset',
                      'XSHOOTER_slt_obs_GenericOffset',
                      'XSHOOTER_slt_obs_Mapping',   
                      'XSHOOTER_slt_cal_SpecphotStdNodding',
                      'XSHOOTER_slt_cal_SpecphotStdOffset',
                      'XSHOOTER_slt_cal_SpecphotStdStare',
                      'XSHOOTER_slt_cal_TelluricStdNod',
                      #'XSHOOTER_slt_cal_TelluricStdOffset', non-existent ???
                      'XSHOOTER_slt_cal_TelluricStdStare',
                      # IFU
                      'XSHOOTER_ifu_acq',
                      'XSHOOTER_ifu_acq_rrm',
                      'XSHOOTER_ifu_obs_GenericOffset',
                      'XSHOOTER_ifu_obs_FixedSkyOffset',
                      'XSHOOTER_slt_obs_Mapping', 
                      'XSHOOTER_ifu_cal_SpecphotStdOffset',
                      'XSHOOTER_ifu_cal_SpecphotStdStare',
                      'XSHOOTER_ifu_cal_TelluricStdOffset',
                      'XSHOOTER_ifu_cal_TelluricStdStare',
                      # IMG
                      'XSHOOTER_img_acq',
                      'XSHOOTER_img_acq_FlatSky',
                      'XSHOOTER_img_obs',
                      'XSHOOTER_img_obs_GenericOffset', 
                      ]

# The acquisition and science parameters that matter for the finding charts.

xshooter_acq_params = {'bos_ra': 0,      # blind offset RA
                       'bos_dec': 0,     # blind offset Dec
                       'acq_pa': 0,      # PA
                       'is_gs': False,   # Guide star defined by user ?
                       'gs': None,       # gs coord as SkyCoord
                  }    
         
xshooter_sci_params = {'ins_mode': None, # Can be 'slt','ifu','img' depending on template
                       'uvb_slt': 0*u.arcsec,
                       'vis_slt': 0*u.arcsec,
                       'nir_slt': 0*u.arcsec,
                       'slt_throw': 0*u.arcsec,
                       'noff': 1,        # Number of offsets
                       'obstype': ['O'], # Observation type (S,O)
                       'coordtype': 'SKY', # Offset type (SKY, DETECTOR) 
                       # NOTE: XSHOOTER allows for SLIT too, but I just call this DETECTOR as well
                       'off1': [0],      # RA offsets
                       'off2': [0],      # DEC offsets
                       'fixoff1': None,  # Fixed RA offsets to sky
                       'fixoff2': None,  # Fixed RA offsets to sky
                       'return': True,   # return to origin ?
                       }
               

# ----------------------------------------------------------------------------------------
def detector_to_sky(dx,dy,pa):
   '''
   Converts XSHOOTER offsets from the DETECTOR reference frame to the SKY reference frame.
   
   Args:
      dx (float): offset in arcsec
      dy (float): offset in arcsec
      pa (float): current position angle in degrees measured East-of-North (XSHOOTER convention)
        
   Returns:
      tuple of floats: the (dra,ddec) offsets in SKY convention.
      
   '''
   
   # Rotate the offsets in the Ra-Dec frame
   # Issue #7: https://github.com/fpavogt/fcmaker/issues/6
   ddec = - np.cos(np.radians(pa)) * dx + np.sin(np.radians(pa)) * dy 
   dra = np.sin(np.radians(pa)) * dx + np.cos(np.radians(pa)) * dy
   
   # Flip the RA axis to increase along East
   dra *= -1.
   
   return (dra, ddec)

# ----------------------------------------------------------------------------------------
def get_p2fcdata_xshooter(fc_params, ob, api):
   '''
   Extracts all the important info to build a finding chart from a given XSHOOTER OB from p2.
   
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
      if not(t['templateName'] in xshooter_templates):
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
         fc_params['acq'] = copy.deepcopy(xshooter_acq_params)
         
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
               
            # Position Angle
            if param['name'] == 'TEL.ROT.OFFANGLE':      
               fc_params['acq']['acq_pa'] = param['value']    
            
         # Store the GS and TTS as SkyCoords
         if fc_params['acq']['is_gs']:
            fc_params['acq']['gs'] = SkyCoord(gs_ra, gs_dec, 
                                       obstime = fcm_m.obsdate, 
                                       #equinox=ob['target']['equinox'], 
                                       frame='icrs',unit=(u.hourangle, u.deg))
      
         # Deal with the blind offset by definition, the "target" in the OB is the 
         # Acq. location, to which we "add" the blind offset to go to the real Target
         fc_params['target'] = fcm_t.offset_coord(fc_params['target'], 
                                                  delta_ra=fc_params['acq']['bos_ra']*u.arcsec, 
                                                  delta_dec=fc_params['acq']['bos_dec']*u.arcsec)
      
         if fc_params['acq']['acq_pa'] == 9999:
            fc_params['tags'] += ['parallactic_angle']
      
            if fcm_m.do_parang:
               # Ok, I need to compute the parallactic angle !
               # Use astroplan to do that
               UT2 = Observer(location=fcm_m.UT2_loc)
               fc_params['acq']['acq_pa'] = UT2.parallactic_angle(fcm_m.obsdate, fc_params['target']).to(u.deg).value 
      
      elif t['type'] in ['science','calib']:
         
         # Ok, I found one more Science template
         fc_params['n_sci'] += 1
         
         temp_name = 'sci%i' %(fc_params['n_sci'])
         
         # Store all the relevant parameters in a dictionary
         fc_params[temp_name] = copy.deepcopy(xshooter_sci_params)
         
         tpl,tplVersion = api.getTemplate(obId, t['templateId'])
         
         # Set my own instrument mode, so I know what I am looking at later on.
         if 'ifu' in t['templateName']:
            fc_params[temp_name]['ins_mode'] = 'ifu'
         elif 'slt' in t['templateName']:
            fc_params[temp_name]['ins_mode'] = 'slt'
         elif 'img' in t['templateName']:
            fc_params[temp_name]['ins_mode'] = 'img'
         else:
            raise Exception('Ouch! This error is impossible!')
            
         for param in tpl['parameters']:

            if param['name'] == 'INS.OPTI3.NAME':             
               fc_params[temp_name]['uvb_slt'] = float(param['value'].split('x')[0])*u.arcsec
            if param['name'] == 'INS.OPTI4.NAME':             
               fc_params[temp_name]['vis_slt'] = float(param['value'].split('x')[0])*u.arcsec
            if param['name'] == 'INS.OPTI5.NAME':             
               fc_params[temp_name]['nir_slt'] = float(param['value'].split('x')[0])*u.arcsec 
            if param['name'] == 'SEQ.NOD.THROW':
               fc_params[temp_name]['slt_throw'] = param['value']*u.arcsec 
            if param['name'] == 'SEQ.NOFFSET':             
               fc_params[temp_name]['noff'] = param['value']
            if param['name'] == 'SEQ.OBS.TYPE':             
               # Warning, here, XSHOOTER differs from MUSE !!!
               fc_params[temp_name]['obstype'] = [i for i in param['value'].split(' ')]
            if param['name'] == 'SEQ.OFFSET.COORDS':
               # NOTE: XSHOOTER allows for SLIT too, but I just call this DETECTOR as well           
               if param['value'] =='SLIT':
                  fc_params[temp_name]['coordtype'] = 'DETECTOR'
               else: 
                  fc_params[temp_name]['coordtype'] = param['value'] 
                
            if param['name'] == 'SEQ.RELOFF1':               
               fc_params[temp_name]['off1'] = param['value']
            if param['name'] == 'SEQ.RELOFF2':             
               fc_params[temp_name]['off2'] = param['value']
            if param['name'] == 'SEQ.OFFSET.ZERO':
               fc_params[temp_name]['return'] = param['value']
            if param['name'] == 'SEQ.FIXOFF.RA':
               fc_params[temp_name]['fixoff1'] = param['value']
            if param['name'] == 'SEQ.FIXOFF.DEC':
               fc_params[temp_name]['fixoff2'] = param['value']
            
         
         # Deal with the AutoNodOnSlit template
         # Pretend the nodding are regular offsets. Ignore any jittering.
         if fc_params[temp_name]['slt_throw'].value >0:
            fc_params[temp_name]['noff'] = 2
            fc_params[temp_name]['obstype'] = ['O','O']
            fc_params[temp_name]['coordtype'] = 'DETECTOR'
            fc_params[temp_name]['off1'] = [0,0]
            fc_params[temp_name]['off2'] = [-fc_params[temp_name]['slt_throw'].to(u.arcsec).value/2.,    
                                            +fc_params[temp_name]['slt_throw'].to(u.arcsec).value]
         
         # Now, deal with the FixedSkyOffset templates.
         # Basically "converts them to a normal GenericOffset, and ignore any jitter
         if not(fc_params[temp_name]['fixoff1'] is None) and not(fc_params[temp_name]['fixoff2'] is None):
            fc_params[temp_name]['noff'] = 2
            fc_params[temp_name]['obstype'] = ['O','S']
            fc_params[temp_name]['coordtype'] = 'SKY'
            fc_params[temp_name]['off1'] = [0,fc_params[temp_name]['fixoff1']]
            fc_params[temp_name]['off2'] = [0,fc_params[temp_name]['fixoff2']]
                      
   return fc_params

# ------------------------------------------------------------------------------
def get_localfcdata_xshooter(fc_params,inpars):
   '''
   Extracts all the important info to build a finding chart from a given XSHOOTER OB defined
   locally.
   
   Args:
      inpars: A dictionnary containing the OB parameters
   Returns:
       A dictionnary containing the MUSE OB parameters
   '''                                    
         
   # Acquisition
   fc_params['acq'] = copy.deepcopy(xshooter_acq_params)
   # Override the user input if clash with instrument mode  
   fc_params['acq']['bos_ra'] = inpars['bos_ra']
   fc_params['acq']['bos_dec'] = inpars['bos_dec']
   fc_params['acq']['is_gs'] = inpars['is_gs']
   fc_params['acq']['gs'] = SkyCoord(inpars['gs_ra'],inpars['gs_dec'], 
                                     frame = 'icrs', 
                                     obstime = Time(fcm_m.obsdate),
                                     equinox = 'J2000')
                                                               
   fc_params['acq']['acq_pa'] = inpars['acq_pa']
   
   
   # Deal with the blind offset by definition, the "target" in the OB is the 
   # Acq. location, to which we "add" the blind offset to go to the real Target
   fc_params['target'] = fcm_t.offset_coord(fc_params['target'], 
                                             delta_ra=fc_params['acq']['bos_ra']*u.arcsec, 
                                             delta_dec=fc_params['acq']['bos_dec']*u.arcsec)
   
   if fc_params['acq']['acq_pa'] == 9999:
      fc_params['tags'] += ['parallactic_angle']
      
      if fcm_m.do_parang:
         # Ok, I need to compute the parallactic angle !
         # Use astroplan to do that
         UT2 = Observer(location=fcm_m.UT2_loc)
         fc_params['acq']['acq_pa'] = UT2.parallactic_angle(fcm_m.obsdate, fc_params['target']).to(u.deg).value 
      
         # Note: if I don't want to show the FoV when parang = 9999, I leave it as is, and 
         # deal with it get_polygon().

   # Observation
   fc_params['n_sci'] = 1
   fc_params['sci1'] = copy.deepcopy(xshooter_sci_params)
   fc_params['sci1']['noff'] = inpars['noff']
   fc_params['sci1']['ins_mode'] = inpars['ins_mode']
   fc_params['sci1']['obstype'] = [i for i in inpars['obstype'][0].split(' ')]
   fc_params['sci1']['off1'] = [float(i) for i in str(inpars['off1'][0]).split(' ')]
   fc_params['sci1']['off2'] = [float(i) for i in str(inpars['off2'][0]).split(' ')]
   fc_params['sci1']['return'] = inpars['return']
   if inpars['coordtype'] == 'SLIT':
      fc_params['sci1']['coordtype'] = 'DETECTOR' # In fcmaker, SLIT = DETECTOR
   else:
      fc_params['sci1']['coordtype'] = inpars['coordtype']

   fc_params['sci1']['uvb_slt'] = inpars['uvb_slt']*u.arcsec
   fc_params['sci1']['vis_slt'] = inpars['vis_slt']*u.arcsec
   fc_params['sci1']['nir_slt'] = inpars['nir_slt']*u.arcsec
   fc_params['sci1']['slt_throw'] = inpars['slt_throw']*u.arcsec
   
   # Deal with the AutoNodOnSlit template
   # Pretend the nodding are regular offsets. Ignore any jittering.
   if fc_params['sci1']['slt_throw'].value >0:
      fc_params['sci1']['noff'] = 2
      fc_params['sci1']['obstype'] = ['O','O']
      fc_params['sci1']['coordtype'] = 'DETECTOR'
      fc_params['sci1']['off1'] = [0,0]
      fc_params['sci1']['off2'] = [-fc_params['sci1']['slt_throw'].to(u.arcsec).value/2.,    
                                            +fc_params['sci1']['slt_throw'].to(u.arcsec).value]

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

   # Define the center of the finding chart to be on the acquisition field.
   target = fc_params['target']  
   
   # Observing fields : build a dictionary for all of them
   fields = {}

   # First, the acquisiton frame 
   counter = 1
   acq_field =  fcm_t.offset_coord(target, 
                                   delta_ra=-fc_params['acq']['bos_ra']*u.arcsec, 
                                   delta_dec=-fc_params['acq']['bos_dec']*u.arcsec)
   # NOTE: I go here from Target to acq, the reverse of the OB convention, hence the '-'
   
   # Start filling the dictionary of fields
   fields[1] = [fc_params['inst'], # Instrument first
                'img', # The instrument mode
                acq_field, # Field central coordinates
                fc_params['acq']['acq_pa'], # Position Angle
                'Acq', # Nature of the field ('Acq', 'Target', 'O', 'S')
                xshooter_acqcam, # The size of the FoV 
                ]
      
   # Include the Target as a "field"
   counter+=1
   fields[counter] = [fc_params['inst'], # Instrument first
                      None,
                      target, # Field central coordinates
                      fc_params['acq']['acq_pa'], # Position Angle
                      'Target', # Nature of the field ('Acq', 'Target', 'O', 'S')
                      None, # The ins_mode
                      ]  

   # Now, loop through all the Science templates I have
   for n in range(fc_params['n_sci']):

      # What is the name of this specific template ?
      temp_name = 'sci%i' %(n+1)
      
      # Ok, what are the offset sequences
      ins_mode = fc_params[temp_name]['ins_mode']
      off1 = fc_params[temp_name]['off1']
      off2 = fc_params[temp_name]['off2']
      obstype = fc_params[temp_name]['obstype']
      coordtype = fc_params[temp_name]['coordtype']
      
      # Get the correct FoV depending on the ins_mode
      if ins_mode == 'img':
         fov = copy.deepcopy(xshooter_acqcam)
      elif ins_mode == 'ifu':
         fov = copy.deepcopy(xshooter_ifu)
      elif ins_mode == 'slt':
         fov = [np.max([fc_params[temp_name]['uvb_slt'].to(u.deg).value,
                        fc_params[temp_name]['vis_slt'].to(u.deg).value,
                        fc_params[temp_name]['nir_slt'].to(u.deg).value])*u.deg,
                xshooter_slt_length]
      else:
         raise Exception('Ouch! ins_mode unknown: %s' % (ins_mode))
      
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
            
         elif coordtype == 'DETECTOR':
            # Transform the dx, dy values in dra, ddec, given the current PA
            (this_dra,this_ddec) = detector_to_sky(off1[o],off2[o],fc_params['acq']['acq_pa'])
            
            # Then sum these and apply any PA shift
            delta_ra += this_dra
            delta_dec += this_ddec
         
         else:
            raise Exception('Ouch! coordtype unknwown: %s' % (coordtype))
           
            
         # Create the field entry
         fields[counter] = [fc_params['inst'],
                            ins_mode,
                            fcm_t.offset_coord(target, 
                                               delta_ra = delta_ra*u.arcsec,
                                               delta_dec = delta_dec*u.arcsec,),
                            fc_params['acq']['acq_pa'], # Position Angle
                            obstype[o],
                            fov,
                           ]
   
      # Ok, I return to the origin after the end of the OB, so let's reset all the offsets
      # to match the acquisition.
      if fc_params[temp_name]['return']:
         delta_ra = 0
         delta_dec = 0
   
   return fields
   
# ------------------------------------------------------------------------------
def get_polygon(central_coord, pa, fov):
   '''
   Given the central location and position of a field, build a polygon to feed to 
   matplotlib down the line.
   
   Args:
      central_coord: an astropy.SkyCoord entry with the center of the XSHOOTER field
      pa: the position angle of the field (p2 convention)
      fov: length 2 list, with the size of the Field-of-View (in angular astropy.units)
   Returns:
      A list of 2-D coordinates for each corner of the field-of-view.
   '''

   # Very well, having done this, I will also build a list of polygons for each field, 
   # that I can feed directly to matplotlib
   polygon = [] #it's a list
   
   # If I got a parallactic angle, then don't show it.
   if pa == 9999:
      return []
   
   xw = fov[0].to(u.deg).value
   yw = fov[1].to(u.deg).value
   
   corners = np.array([[-xw/2,-yw/2],[-xw/2,yw/2],[xw/2,yw/2],[xw/2,-yw/2]])                       
  
   # Rotate the different corners by the correct amount
   for (j,corner) in enumerate(corners):    
      corners[j] = np.dot(corners[j],fcm_t.myrotmatrix(pa))
    
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
   The specific XSHOOTER function that draws a specific observation field.
   
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
                 'lw':0.75, 'ms':50, 'ls':':', 'mc':'royalblue'},
            'S':{'c':'darkcyan', 'lwm':1, 'zorder':5, 'marker':'s', 
                 'lw':0.75, 'ms':50, 'ls':':', 'mc':'darkcyan'},
            }
   
   
   this_coords = [field[2].ra, field[2].dec]
   
   for ax in [ax1,ax2]:
      # Show the center of the field, without obstructing it to see the blind offset star
      if not(field[4] in['Acq','Target']) or (ax is ax1):
         ax.show_markers(this_coords[0].deg, this_coords[1].deg, 
                         marker=skins[field[4]]['marker'],
                         edgecolor=skins[field[4]]['mc'],
                         s=skins[field[4]]['ms'], 
                         linewidth=skins[field[4]]['lwm'], 
                         zorder=skins[field[4]]['zorder'],
                         ) 
      
      # instrument footprint, except for the "target"
      if field[4] != 'Target':
         ax.show_polygons(get_polygon(field[2],field[3],field[5],), 
                           edgecolor = skins[field[4]]['c'],
                           linewidth = skins[field[4]]['lw'],
                           zorder = skins[field[4]]['zorder'],
                           linestyle = skins[field[4]]['ls'],)
           
      # the GS validity area
      ax.show_circles([this_coords[0].deg],
                      [this_coords[1].deg],
                      [((fcm_m.outer_GS_Cas*u.arcsec).to(u.degree)).value,],
                       color='k', lw=0.5, 
                       zorder=skins[field[4]]['zorder']
                     )           
                                    
                   
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
            if ((field[2].separation(fc_params['acq']['gs']) > (fcm_m.outer_GS_Cas*u.arcsec)) or \
                (field[2].separation(fc_params['acq']['gs']) < (inner_GS_search*u.arcsec))):
                
               if fcm_m.fcm_usetex:
                  lab = r'\textbf{!}'
               else:
                  lab = '!'
                  
               ax.add_label(fc_params['acq']['gs'].ra.deg,fc_params['acq']['gs'].dec.deg,lab, 
                            verticalalignment='center', 
                            horizontalalignment='center',size=12,color='crimson',
                            bbox=dict(boxstyle="circle,pad=0.17", facecolor='w',ec='crimson', alpha=1))
       
      # Here, show better the orientation of the slit
      if (field[1] == 'slt') and (field[3] != 9999):
         #ax1.add_label(0.9,0.9, '--', relative=True, horizontalalignment='center', size=70,
         #              verticalalignment='center', 
         #              rotation = 90+field[3], color=skins['Acq']['mc'])
         ax1.add_label(0.9,0.9, '--- p.a. ---', relative=True, horizontalalignment='center', size=10,
                       verticalalignment='center', 
                       rotation = 90+field[3], color=skins['Acq']['mc'],
                       bbox=dict(boxstyle="round4,pad=0.05", facecolor='w',ec='w', alpha=1))
                       
      
          
       
          
   
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
                    
                    
                    
                    