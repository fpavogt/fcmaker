# -*- coding: utf-8 -*-
'''
 fcmaker: a Python module to automatically create finding charts for ESO OBs in p2.\n
 Copyright (C) 2017,  F.P.A. Vogt
 --- oOo ---
 
 This file contains the master fcmaker routines. See the dedicated website for more info:
 http://fpavogt.github.io/fcmaker
 Created October 2017, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
 
 --- oOo ---
  
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
'''
# ----------------------------------------------------------------------------------------
# Import the fcmaker submodules and things -----------------------------------------------
from . import fcmaker_metadata as fcm_m
from .fcmaker_metadata import __version__
from . import fcmaker_tools as fcm_t
from . import fcmaker_plots as fcm_p
from . import fcmaker_instrument_dispatch as fcm_id
from . import fcmaker as fcm

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
from . import p2api

# Use argparse to make fcmaker user friendly ---------------------------------------------
parser = argparse.ArgumentParser(description='''Create ESO-compliant finding charts from \
                                                p2 (or locally) from a parameter file. \
                                                If no parameter file is passed, connect \
                                                to p2 and select an obId manually. ''',
                                 epilog =' Full documentation: %s \n \n \
                                           Feedback, questions, comments: \
                                           frederic.vogt@alumni.anu.edu.au \n' % (os.path.join(fcm_m.fcm_dir,'docs','index.html')))

parser.add_argument('--version', action='version', version=('fcmaker %s'%__version__))

parser.add_argument('-l','--local', action='store_true',
                    help='feed a manual, local OB description')

parser.add_argument('--clear-SkyView-cache', action='store_true',
                    help='clear the SkyView cache')
                                 
parser.add_argument('--no-montage', action='store_true',
                    help='disable the use of Montage to rotate the charts')   
                                                  
parser.add_argument('--no-systemtex', action='store_true',
                    help='disable the use of the system-wide LaTeX')

parser.add_argument('--no-upload', action='store_true',
                    help='disable the uload of finding chart to p2')                    

parser.add_argument('--do-pdf', action='store_true',
                    help='save a pdf version of the chart (in addition to the jpg)')                    

parser.add_argument('--do-png', action='store_true',
                    help='save a png version of the chart (in addition to the jpg)') 

parser.add_argument('--obsdate', action='store', metavar='obsdate', nargs='+', 
                    default=[datetime.strftime(datetime.now(), '%Y-%m-%d')],
                    help='Date of the observations (for targets with proper motions)')                       
                                                             
# An input parameter file - force to use only 1 file at a time                                      
parser.add_argument('-f', action='store', metavar='filename', nargs=1, 
                    type=open, help='parameter filename')
                    

# Start of the interactive part ----------------------------------------------------------
if __name__ == "__main__":

   # What did the user type in ?
   args = parser.parse_args()
   
   # Open the docs ?
   '''
   if args.show_doc_loc:
   
      # locate the docs
      fn = os.path.join(fcm_m.fcm_dir,'docs','index.html')
      if not(os.path.isfile(fn)):
         raise Exception('Docs not found. %s is not a file!' %(fn))
      else:
         print('Doc location: %s' %(fn)) 
   '''
   
   # Avoid using the system-wide latex for the plots ? -----------------------------------
   fcm_m.fcm_usetex = not(args.no_systemtex)
   if args.no_systemtex:  
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
   
   # Set some of the other variables -----------------------------------------------------
   # Use Montage to force the plots North?
   fcm_m.set_North = not(args.no_montage)
      
   # Clear SkyView Cache ?
   fcm_m.clear_SkyView_cache = args.clear_SkyView_cache
   
   # What is the proposed time of the observation ?
   # Either specified by the User, or default to the current time.
   fcm_m.obsdate = dup.parse(' '.join(args.obsdate), default=datetime(2018, 7, 1, 0, 0))
   
   # If not specified, assume UTC time zone   
   if fcm_m.obsdate.tzinfo is None:
      #warnings.warn(' "--obsdate" timezone not specified. Assuming UTC.')
      fcm_m.obsdate = fcm_m.obsdate.replace(tzinfo=pytz.utc)
         
   # Did I get a parameter file ? --------------------------------------------------------
   if not(args.f is None): # Yes !
      
      # Load the parameter file
      inpars = yaml.load(args.f[0])
      
      # If I need to connect to p2, extract the user ID and obIDs. -----------------------
      # If no password or user ID in file, ask for it now.
      # If local file, assign "fake obId"
      if args.local:
         obIDs = [1]
      
      else:
         userID = inpars['userId']
         obIDs = inpars['obIds']
      
         if inpars['userId'] is None:
            userID = input('p2 user Id: ')
         else:
            userID = inpars['userId']

         if inpars['pswd'] is None:
            pswd = getpass.getpass('Password: ')
         else:
            pswd = inpars['pswd']
      
      # Set the relative paths to store the data as requested ----------------------------
      if 'data_loc' in inpars.keys():
         fcm_m.data_loc = os.path.join('.',inpars['data_loc'])
      if 'plot_loc' in inpars.keys():
         fcm_m.plot_loc = os.path.join('.',inpars['plot_loc']) 
         
      # Did the user specify custom background images ? ----------------------------------
      if 'bk_images' in inpars.keys():
         
         if args.local:
            bk_images = [inpars['bk_images']]
         
         else:
            # Make sure I have a bk_image specified for each ob Id !
            if (type(inpars['bk_images']) == list) and (len(inpars['bk_images']) != len(obIDs)):
               raise Exception('Ouch! Please specify one bk_image for each obID!')
               
            elif (type(inpars['bk_images']) == list): # A different background image has been specified for each OB
               bk_images = inpars['bk_images']
               
            else: # Ok, we want the same background image for all ... let's duplicate it
               bk_images = []
               for i in range(len(obIDs)):
                  bk_images += [inpars['bk_images']]
               
               
      # Idem for the associated wavelengths     
      if 'bk_lams' in inpars.keys():
         
         if args.local:
            if not(inpars['bk_lams'] == 'None'):
               bk_lams = [inpars['bk_lams']] 
            else:
               bk_lams = [None]
         
         else:
            # Make sure I have a bk_lam specified for each ob Id !
            if (type(inpars['bk_lams']) == list) and (len(bk_lams) != len(obIDs)):
               raise Exception('Ouch! Please specify one bk_lam for each obID!')
            
            elif (type(inpars['bk_lams']) == list) :
                bk_lams = inpars['bk_lams']
            
            else:
               bk_lams = []
               for i in range(len(obIDs)):
                  if not(inpars['bk_lams']== 'None'):
                     bk_lams += [inpars['bk_lams']]
                  else:
                     bk_lams += [None]        
            
      
   else: # ok, a manual run it is ! ------------------------------------------------------

      # Enforce the use of parameter file if creating a local finding chart 
      if args.local:
         raise Exception('Ouch! To use the local mode, you must specify a suitable '+
                          'parameter file with "-f filename".')

      # Get the p2 login info
      userID = input('Enter your p2 userID: ')
      pswd = getpass.getpass('Password: ')

      # What OB should we create a finding hart for ?
      obIDs = [eval(input('OB Id from p2: '))]
      
      # Just use the default bk_image in manual mode
      bk_images = [None]
      bk_lams = [None]
  
    
   # Make sure I have the necessary folders to store stuff -------------------------------
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

   #  Log into p2 if required  -----------------------------------------------------------
   if args.local:
      api = None
   else:
      api = p2api.ApiConnection('production',userID,pswd, False) 
      
      # Clear the password right away, for security reasons
      pswd = None
      
   #  Now, start doing stuff ----------------------------------------------------------------                      
   for (o,obID) in enumerate(obIDs):
      
      print(' ') # Some space for clarity
      
      # Step 1: extract the OB parameters ------------------------------------------------
      if args.local:
         fc_params = fcm_id.get_localfcdata(inpars)#local_OB
      else:
         print('%i: fetching the OB parameters ...' % (obID))
         # Fetch all the parameters required for the finding chart
         fc_params = fcm_id.get_p2fcdata(obID, api)
   
      # Step 2: create the finding chart -------------------------------------------------
      print('%i: creating the finding chart ...' % (obID))
      # Send these to the plotting routine
      fc_fn = fcm_p.make_fc(fc_params, bk_image=bk_images[o], bk_lam=bk_lams[o],
                            do_pdf=args.do_pdf, do_png=args.do_png)
   
      # If I am working locally, I am done.
      if args.local:
         continue
   
      # Step 3: upload the chart to p2 ---------------------------------------------------
      # If the no_upload flag is set, keep going right away
      if args.no_upload:
         continue
      
      answer = None
      while not(answer in ['y','n']):
         answer = input('%i: attach finding chart to p2 OB (y/n)? ' % (obID))
      
      # Ok, let's start looking at the finding charts online 
      if answer == 'y':
         
         # List existing finding charts already attached to the OB
         fcNames, _ = api.getFindingChartNames(obID)
         
         print('  Existing finding charts:')
         for i in range(5):
            if i < len(fcNames):
               print('  %i: %s' % ((i+1),fcNames[i]) )
            else:
               print('  %i: empty' % ((i+1)) )
         answer = None
         while not(answer in range(0,6)):
            answer = eval(input('  Which slot to upload to (1-5; 0 = no upload)? '))
         
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
            print('  Done.')
         else:
            print('  Ok, moving on ...')
            
   
# ----------------- End of the World as we know it ---------------------------------------
# ----------------------------------------------------------------------------------------