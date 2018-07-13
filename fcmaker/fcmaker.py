# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------------------

import numpy as np
import warnings

from astropy.coordinates.sky_coordinate import SkyCoord
from astropy import units as u

# Import the sub-routines
from . import fcmaker_metadata as fcm_m
from . import fcmaker_instrument_dispatch as fcm_id
from .fcmaker_metadata import __version__
from . import fcmaker_tools as fcm_t
from . import fcmaker_plots as fcm_p

# Import generic python packages
import os
import sys
import numpy as np
import getpass
import yaml
import argparse
import webbrowser
from datetime import datetime
from dateutil import parser as dup
import pytz
import warnings

# Some astrostuff
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy import units as u

# Import matplotlib to set the proper style file here already
from matplotlib import pylab as plt

# Import some p2 magic 
import p2api



#from . import fcmaker_muse as fcm_muse

'''
 fcmaker: a Python module to automatically create finding charts for ESO OBs in p2.\n
 Copyright (C) 2017,  F.P.A. Vogt
 --- oOo ---
 This file contains general functions for the fcmaker routines. 
 Created October 2017, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''

# ----------------------------------------------------------------------------------------
def set_obsdate(obsdate):
   '''
   A function that sets the obsdate parameter inside fcmaker.
   
   Args:
      obsdate: string. The observing time (and date) for the OB. Should be parsable by dateutils.parser.parse(). None = now.
   '''
   
   # What is the proposed time of the observation ?
   if obsdate is None:
      obsdate = datetime.strftime(datetime.utcnow(), '%Y-%m-%d %H:%M:%S') + ' UTC'
   
   # Store it
   fcm_m.obsdate = dup.parse(obsdate)
   
   # If not specified, assume UTC time zone   
   if fcm_m.obsdate.tzinfo is None:
      warnings.warn(' [fcmaker] obsdate timezone not specified. Assuming UTC.')
      fcm_m.obsdate = fcm_m.obsdate.replace(tzinfo=pytz.utc)
         
# ----------------------------------------------------------------------------------------
def set_systemtex(systemtex):
   '''
   A function that sets whether fcmaker uses the default system LateX, or not. 
   
   Args:
      no_systemtex: bool. True = use Python LaTex. False = System LaTeX.
   '''
   
   # Avoid using the system-wide latex for the plots ?
   fcm_m.fcm_usetex = systemtex
   if not(systemtex):  
      fcm_m.fcm_plotstyle = os.path.join(fcm_m.fcm_dir,'mpl_styles',
                                         'fcmaker_plots_nolatex.mplstyle')
      plt.style.use(fcm_m.fcm_plotstyle)
      # I don't like this next line too much ... but failed to find a better way!
      #fcm_p.reload() 
   else:
      fcm_m.fcm_plotstyle = os.path.join(fcm_m.fcm_dir,'mpl_styles',
                                         'fcmaker_plots.mplstyle')
      plt.style.use(fcm_m.fcm_plotstyle)
      # I don't like this next line too much ... but failed to find a better way!
      #fcm_p.reload() # I don't like this reload too much either ...

# Function to setup fcmaker parameters ---------------------------------------------------   
def set_fcmaker(systemtex,montage,clear_SkyView_cache, data_loc, plot_loc, do_parang):
   '''
   A function that sets the generic options of fcmaker. 
   
   Args:
      no_systemtex: bool. True = use Python LaTex. False = System LaTeX.
      no_montage: bool. True = do NOT use Montage.
      clear_SkyView_cache: bool. True = clear the cache.
   '''
   # Set system LateX or not
   set_systemtex(systemtex)
   
   # Use Montage to force the plots North?
   fcm_m.set_North = montage
      
   # Clear SkyView Cache ?
   fcm_m.clear_SkyView_cache = clear_SkyView_cache

   # Set the relative paths to store the data as requested
   fcm_m.data_loc = os.path.join('.',data_loc)
   fcm_m.plot_loc = os.path.join('.',plot_loc) 
   
   # Make sure I have the necessary folders to store stuff
   # Create them if necessary.
   for folder in [fcm_m.data_loc, fcm_m.plot_loc]:
      if not(os.path.isdir(folder)):
         answer = None
         while not(answer in ['y','n']):
            answer = input('Creating local directory ./%s (y/n)? ' % (folder))
            
         if answer =='y':
            os.mkdir(folder)
         else:
            raise Exception('Ouch! I really wanted to create that folder!')

   # Set whether I want to deal with parallactic angles (or not)
   fcm_m.do_parang = do_parang

# Function to create FC from p2 ----------------------------------------------------------
def make_fc( p2uid = None, pswd = None,  
             obids = [],
             bk_images = [],
             bk_lams = [],     
             data_loc = os.path.join('.','fcm_data'),        
             plot_loc = os.path.join('.','fcm_plots'),
             do_pdf = True, do_png = False,
             no_upload = False, 
             systemtex = False,
             montage = False,
             clear_SkyView_cache = False,  
             obsdate = None,
             do_parang = False,
             ):
   '''
   The main fcmaker function, to create finding charts from p2.
   
   Args:
      p2uid: string. P2 user ID (will prompt if is None)
      pswd:  string. P2 user password (will prompt if is None)
      obids: list of int. List of the P2 OB ID for which to generate finding charts
      bk_images: list of string. Specify the FC background: SkyView survey name, None for default, or local FITS filename
      bk_lams: list of string. Specify the wavelength of the background charts. None to read from FITS header.
      data_loc: string. Relative path to the background images (local FITS files or SkyView images).
      plot_loc: string. Relative path to the background images (local FITS files or SkyView images).
      do_pdf: bool. Save a local PDF file ?
      do_png: bool. Save a local png file ?
      no_upload: bool. Skip the upload of finding charts to p2 ?
      systemtex: bool. Use the system Latex instead of the Python latex ?
      montage: bool. Use of Montage to rotate the fields with North up ?
      clear_SkyView_cache: bool. Clear the SkyView cache ?
      obsdate: string. Year (Month, Day, Hour, Minute, ...) of the observation
      do_parang: bool. Show the instrument field-of-view when a parallactic angle is required ?
   '''
 
   starttime = datetime.now()
   
   # Make sure I have at least one obids listed if not in local mode 
   if len(obids) ==0:
      obids = [eval(input('OB Id from p2: '))]

   # Make some safety checks 
   if not(type(bk_images) == list):
      raise Exception('Ouch! bk_images must be a list !')
   if not(type(bk_lams) == list):
      raise Exception('Ouch! bk_lams must be a list !')
   if not(type(obids) == list):
      raise Exception('Ouch! obids must be a list !')
   
   # Set the observing date (and time) 
   set_obsdate(obsdate) 
   
   # Set the generic parameters for fcmaker
   set_fcmaker(systemtex,montage,clear_SkyView_cache, data_loc, plot_loc, do_parang) 
  
   # If I need to connect to p2, extract the user ID and obIDs. 
   # If no password or user ID in file, ask for it now.
   if p2uid is None:
         p2uid = input('p2 user ID: ')

   if pswd is None:
         pswd = getpass.getpass('Password: ')
         
   # Did the user specify custom background images ? 
   if len(bk_images) >0:
         
      # Make sure I have a bk_image specified for each ob Id !
      if len(bk_images) != len(obids):
         raise Exception('Ouch! Please specify one bk_image for each obID!')
               
   else: 
      bk_images = [None]*len(obids)
                    
   # Idem for the associated wavelengths 
   if len(bk_lams) > 0:    
   
      # Make sure I have a bk_lam specified for each ob Id !
      if len(bk_lams) != len(obids):
         raise Exception('Ouch! Please specify one bk_lam for each obID!')
   
   else:
      bk_lams = [None]*len(obids)    
   
   #  Log into p2
   api = p2api.ApiConnection('production',p2uid,pswd, False) 
      
   # Clear the password right away, for security reasons
   pswd = None
      
   #  Now, start doing stuff                    
   for (o,obID) in enumerate(obids):
      
      print(' ') # Some space for clarity
      
      # Step 1: extract the OB parameters
      print('%i: fetching the OB parameters ...' % (obID))
      # Fetch all the parameters required for the finding chart
      fc_params = fcm_id.get_p2fcdata(obID, api)
   
      # Step 2: create the finding chart 
      print('%i: creating the finding chart ...' % (obID))
      # Send these to the plotting routine
      fc_fn = fcm_p.draw_fc(fc_params, bk_image=bk_images[o], bk_lam=bk_lams[o],
                            do_pdf=do_pdf, do_png=do_png)
   
      # Step 3: upload the chart to p2
      # If the no_upload flag is set, keep going right away
      if no_upload:
         continue
      
      # Ok, let's start looking at the finding charts online 
      # List existing finding charts already attached to the OB
      fcNames, _ = api.getFindingChartNames(obID)
      
      # Are there any finding charts over there ?
      if len(fcNames)>0:
      
         print('  Existing finding charts:')
         for i in range(5):
            if i < len(fcNames):
               print('  %i: %s' % ((i+1),fcNames[i]) )
            else:
               print('  %i: empty' % ((i+1)) )
         answer = None
         while not(answer in range(0,6)):
            answer = eval(input('  Which slot to upload to (1-5; 0 = no upload)? '))
            
      else:
         # If no finding charts exist, then just put it in the first spot (no need to ask)
         print('  No finding chart in the OB (yet): using slot 1 for upload ...')
         answer = 1
               
      # Check if the finding chart slot is busy ... do we want to replace it ?
      if answer == 0:
         fill = 'n'
      elif answer <= len(fcNames):
         fill = None
         while not(fill in ['y','n']):
            fill = input('  Finding chart slot occupied: overwrite (y/n)? ')   
      else:
         fill = 'y'
         
      if fill == 'y':
         # If I have to delete the existing finding chart
         if answer <= len(fcNames):
            api.deleteFindingChart(obID, answer) # p2 finding chart index start at 1
            
         # And upload the new chart
         api.addFindingChart(obID, fc_fn)
   
   print('All finding charts done in %.1f seconds.' %((datetime.now()-starttime).total_seconds()))
   
# Function to create FC from a local file 
def make_fc_local(f, 
                  do_pdf = True, do_png = False,
                  systemtex = False,
                  montage = True,
                  clear_SkyView_cache = False,  
                  obsdate = None,
                  do_parang = False,
                  ):
   '''
   The other fcmaker function, to create finding charts from a local file.
   
   Args:
      f: an open file, e.g. _io.TextIOWrapper from open(filename)
      do_pdf: bool. Save a local PDF file ?
      do_png: bool. Save a local png file ?
      systemtex: bool. Use the system Latex instead of the Python latex ?
      montage: bool. Use of Montage to rotate the fields with North up ?
      clear_SkyView_cache: bool. Clear the SkyView cache ?
      obsdate: string. Year (Month, Day, Hour, Minute, ...) of the observation
      do_parang: bool. Show the instrument field-of-view when a parallactic angle is required ?
   '''
   
   starttime = datetime.now()
   
   if not(os.path.isfile(f.name)):
      raise Exception('Ouch! unknown file: %s' % (fn))
   
   # Load the parameter file
   inpars = yaml.load(f)
   
   # Set the observing date (and time)
   set_obsdate(obsdate)
   
   # Set the generic parameters for fcmaker 
   set_fcmaker(systemtex, montage, clear_SkyView_cache, 
               inpars['data_loc'], inpars['plot_loc'], do_parang)
   
   # Build the fidning chart dictionnary
   fc_params = fcm_id.get_localfcdata(inpars) 
   
   if fc_params['ob_name'] is None:
      raise Exception('Ouch! Please specify an "obname" in the local file')
   
   if fc_params['pi'] is None:
      fc_params['pi'] = ' '
   
   if fc_params['prog_id'] is None:
      fc_params['prog_id'] = ' '
   
   if fc_params['ob_id'] is None:
      fc_params['ob_id'] = -2
   
 
   # Step 2: create the finding chart
   print('Creating the finding chart from local info ...')
   # Send these to the plotting routine
   fc_fn = fcm_p.draw_fc(fc_params, 
                         bk_image=inpars['bk_image'], bk_lam=inpars['bk_lam'],
                         do_pdf=do_pdf, do_png=do_png)
   
   
   print('All done in %.1f seconds.' %((datetime.now()-starttime).total_seconds()))