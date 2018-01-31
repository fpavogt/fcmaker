# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------------------

import os
import numpy as np

import warnings # for aplpy displaying nasty packages
from datetime import datetime

# Import plotting packages
from matplotlib import pylab as plt
import matplotlib.gridspec as gridspec

with warnings.catch_warnings():
   warnings.simplefilter("ignore")
   import aplpy

# Some astro stuff
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy import units as u
import astropy.io.fits as fits
from astropy.constants import c
from astroquery.skyview import SkyView
from astroquery.vizier import Vizier

# My own mess
from . import fcmaker_metadata as fcm_m
from . import fcmaker_tools as fcm_t
from . import fcmaker_instrument_dispatch as fcm_id

# Set the fcmaker plot style
plt.style.use(fcm_m.fcm_plotstyle)
#plt.style.reload_library()

'''
fcmaker: a Python module to automatically create finding charts for ESO OBs in p2.\n
Copyright (C) 2017,  F.P.A. Vogt
--- oOo ---
This file contains general tools for the fcmaker plotting routines.
Created October 2017, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''

# ----------------------------------------------------------------------------------------
def add_orient(ax, field_center, 
               radius = 540*u.arcsec, arrow_width = 120*u.arcsec, usetex=True):
   '''
   Plot some N-E orientation arrows
   
   Rely on the show_arrow function from aplpy, which should be able to deal with
   any image orientation, even is North /= up.
   
   Args:
      ax: the aplpy ax instance to add the arrows to
      field_center: a SkyCoord entry marking the central pointing of the chart
      radius: the distance from the field_center to draw the arrows at, in arcsec.
      arrow_width: the arrow size, in arcsec
      usetex: whether a system LaTeX is used, or not
   
   '''
   # Center of the cross
   aloc = np.array([field_center.ra.deg-radius.to(u.deg).value/np.cos(field_center.dec.radian),
                    field_center.dec.deg-radius.to(u.deg).value])
                    
   a_south = aloc - np.array([0,arrow_width.to(u.deg).value/2.])
   a_west = aloc - np.array([arrow_width.to(u.deg).value/2./np.cos(field_center.dec.radian),0])
   
   ax.show_arrows([a_west[0], a_south[0]], [a_west[1],a_south[1]], 
                  [arrow_width.to(u.deg).value/np.cos(np.radians(aloc[1])),0], 
                  [0,arrow_width.to(u.deg).value], 
                  head_length=20,head_width=15, color='k')
   
   if usetex:
      lab1 = r'\smaller\textbf{N}'
      lab2 = r'\smaller\textbf{E}'
   else:
      lab1 = r'N'
      lab2 = r'E'
         
   ax.add_label(a_south[0],a_south[1]+1.3*arrow_width.to(u.deg).value,lab1,
                    verticalalignment='center',horizontalalignment='center',
                    color = 'k', weight='heavy', size='large',
                    bbox=dict(facecolor='none',ec='none', alpha=0.75, boxstyle="round,pad=.3"))
   ax.add_label(a_west[0]+1.3*arrow_width.to(u.deg).value/np.cos(np.radians(aloc[1])),a_west[1],lab2,
                    verticalalignment='center',horizontalalignment='center',
                    color = 'k', weight = 'heavy', size='large',
                    bbox=dict(facecolor='none',ec='none', alpha=0.75, boxstyle="round,pad=.3"))


# ----------------------------------------------------------------------------------------
def get_bk_image(bk_image, fc_params):
   ''' 
   A function that fetches the background image from the web, or locally if specified
   by the user.
   
   Args:
      bk_image: (string) defines whethe to use a local image (FITS filename) or a SkyView 
                image (survey name). 
      fc_params: the parameters of the finding chart (a dictionnary)
   Returns:
      The absolute path of the background image (possibly after download), and the survey name.
   
   .. note::
     
      For local images, the code excpet them to be stored in `fcm_m.data_loc` folger.
      
   '''
   
      # --------------------------------- BKGD ---------------------------------
   if not(bk_image in [None,'None']):
      # Has the user provided his/her own FITS file ?  
      if not(bk_image in SkyView.survey_dict['overlay_blue']):
         # Very well, I guess I have a custom FITS file.
         fn_bk_image = os.path.join(fcm_m.data_loc, bk_image)
         survey = bk_image.split('.')[0] # assume a filename.fits ...
      else:
         fn_bk_image = ''
         survey = bk_image
   else:
      # Nothing specified ... use the default ...
      fn_bk_image = ''
      survey = fcm_id.get_bk_image(fc_params)
   
   # Clear the Skyview cache ?
   if fcm_m.clear_SkyView_cache:
      print('Clearing SkyView cache...')
      cache_loc = SkyView.cache_location
      items  = os.listdir(cache_loc)
      
      for item in items:
         os.remove(os.path.join(cache_loc,item))
   
   # Do I already have an image ?
   if fcm_m.clear_SkyView_cache or not(os.path.isfile(fn_bk_image)):
   
      # Specify FITS filename
      fn_bk_image = os.path.join(fcm_m.data_loc, fc_params['ob_name'].replace(' ','_')+'_' + 
                                 fcm_id.get_bk_image(fc_params).replace(' ','-') + 
                                 '.fits')
  
      # Use SkyView to get the image 
      print('   Downloading the %s background image with SkyView ...' % (survey))
      path = SkyView.get_images(position=fc_params['target'], survey=[survey], 
                                radius=fcm_m.bk_radius,
                                pixels='%i' %((fcm_m.bk_radius/fcm_m.bk_pix).to(u.arcsec/u.arcsec).value),
                                )
      if len(path)==0:
         raise Exception('No image downloaded via SkyView.get_images!')
      
      # Save it locally for later
      outfits = fits.HDUList([path[0][0]])
      outfits.writeto(fn_bk_image, output_verify='fix', overwrite=True)

     
   return (fn_bk_image, survey)
   
   # --------------------------------- BKGD ---------------------------------

def get_bk_image_lam(fn_bk_image, fc_params):
   '''
   Reads the wavelength of the background image from the header. Assumes SkyView image.
   
   Args:
      fn_bk_image: path to background image FITS file.
       fc_params: the parameters of the finding chart (a dictionnary)
   Returns:
      The image wavelength range as a string.
   '''
   
   # Find the wavelength of the data
   hdu= fits.open(fn_bk_image)
   header = hdu[0].header
   hdu.close()
   
   # Extract the frequency range thanks to SkyView
   freqs = [l.split(' ')[-2].split('-') for l in header['COMMENT'] if 'Bandpass' in l][0]
   # WARNING HERE: I assume the bandpass is in THz in all cases !!!!!
   freqs = [np.int(c.value/(np.int(freq)*1.e12)*1.e9) for freq in freqs]
   bk_lam = r'%i - %i nm (%s)' % (freqs[1],freqs[0],fcm_id.get_bk_image(fc_params))

   return bk_lam

# ----------------------------------------------------------------------------------------
def make_fc(fc_params, bk_image = None, bk_lam = None, do_pdf = False, do_png = False):
   ''' 
   The finding chart master plotting function.
   
   Args:
      fc_params: a dictionnary containing the OB parameters
      bk_image: the background image as string, either a SkyView entry, or local filename
      bk_lam: the wavelength of the chart as a string
      do_pdf (bool): save a pdf file, in addition to jpg
      do_png (bool): save a png file, in addition to jpg
   Returns:
      The finding chart filename as a string
   '''
   
   # Load all the fields dictionaries
   fields = fcm_id.get_fields_dict(fc_params)
   
   # Get the background image in place
   (fn_bk_image, survey) = get_bk_image(bk_image, fc_params)
   
   if bk_lam is None:
      bk_lam = get_bk_image_lam(fn_bk_image, fc_params)
   
   # Start the plotting
   plt.close(1)
   fig1 = plt.figure(1, figsize=(14.17,7)) #14.17 inches, at 50% = 1 full page plot in A&A. 

   ax1 = aplpy.FITSFigure(fn_bk_image, figure=fig1, north=fcm_m.set_North,
                          subplot=[0.12,0.12,0.35,0.72])
   ax1.show_grayscale(invert = True, stretch='linear', pmin = fcm_id.get_pmin(survey))

   ax2 = aplpy.FITSFigure(fn_bk_image, figure=fig1, north=fcm_m.set_North, 
                          subplot=[0.59,0.12,0.38,0.8])
   ax2.show_grayscale(invert = True, stretch='linear', pmin = fcm_id.get_pmin(survey))
 
   # Get the radius and center of the charts
   (left_radius, right_radius) = fcm_id.get_chart_radius(fc_params)
   (left_center, right_center) = fcm_id.get_chart_center(fc_params)
   
   ax1.recenter(left_center.ra,left_center.dec,radius=left_radius/3600.)
   ax2.recenter(right_center.ra,right_center.dec,radius=right_radius/3600.) 
   
   # Query UCAC2 via Vizier over the finding chart area, to look for suitable Guide Stars
   
   print('   Querying UCAC2 via Vizier to look for possible Guide Stars ...')
   Vizier.ROW_LIMIT = 10000
   gs_table = Vizier.query_region(right_center, radius=right_radius*u.arcsec,
                                  inner_radius = fcm_id.get_inner_GS_search(fc_params)*u.arcsec,
                                  catalog="UCAC2")
   gs_table = gs_table[gs_table.keys()[0]]
   
   # Turn the table into a list of SkyCoord, that I can clean as I go.
   # Only select guide stars in the nominal GS mag range
   gs_list = [fcm_t.propagate_pm(SkyCoord(ra=line['RAJ2000'],
                                          dec=line['DEJ2000'],
                                          obstime = fcm_m.obsdate, 
                                          equinox = 'J2000', 
                                          frame = 'icrs',
                                          unit=(u.deg, u.deg)), 
                                 2000., line['pmRA']/1000., line['pmDE']/1000.,) 
              for line in gs_table if fcm_m.gs_mag[0] <= line['UCmag'] <= fcm_m.gs_mag[1]]
   
   # Show the observation footprint
   for f in fields:
      
      # Call the function that will plot all the important stuff for this field.
      fcm_id.plot_field(ax1, ax2, fc_params, fields[f])
      
      # Keep all the possible guide stars in the area.
      gs_list = [star for star in gs_list if 
                 (fields[f][2].separation(star) > (fcm_id.get_inner_GS_search(fc_params)/3600.*u.deg)) 
                 and (fields[f][2].separation(star) < (10./60.*u.deg))]         

   # Show all the suitable Guide Star in the area
   if len(gs_list)>0:
      ax2.show_markers([np.array(star.ra) for star in gs_list],
                       [np.array(star.dec) for star in gs_list], 
                       marker='o', 
                       edgecolor='crimson',
                       s=50, 
                       linewidth=1.0)
   else:
      warnings.warn('Watch out ... no suitable Guide Star found in UCAC2 !')
                
   # Add orientation arrows to the large view plot
   add_orient(ax2, right_center, radius = (right_radius*0.8)*u.arcsec, usetex = fcm_m.fcm_usetex)
    
   # Add a scale bar
   (scl,scll) = fcm_id.get_scalebar(fc_params['inst'])
   
   ax1.add_scalebar(scl) # Length in degrees
   ax1.scalebar.show(scl, label=scll, corner = 'bottom left', 
                     color = 'k', frame=1)
   ax1.scalebar.set_linewidth(2)

   # Fine tune things a bit further, just because I can ...
   for ax in [ax1,ax2]:
      ax.tick_labels.set_xformat('hh:mm:ss')
      ax.axis_labels.set_xpad(10)
      ax.ticks.set_linewidth(1.5)
      ax.set_tick_size(10)
      ax.set_tick_color('k')
      ax.axis_labels.set_xtext('R.A. [J2000]')

   ax1.axis_labels.set_ytext('Dec. [J2000]')
   ax1.axis_labels.set_ypad(-10)
   ax2.axis_labels.set_ytext('')

   # Add the required OB information to comply with ESO requirements ...
   # ... and make the life of the night astronomer a lot easier !
   ax1.add_label(0.0,1.18, 'Run ID: '+fc_params['progId']+' | '+fc_params['pi'], 
                 relative=True, horizontalalignment='left', size=16)
   
   # Fix some bugs for anyone not using LaTeX ... sigh ...
   if fcm_m.fcm_usetex:
      lab = fc_params['ob_name'].replace('_','\_')
   else:
      lab = fc_params['ob_name']             
                
   ax1.add_label(0.0,1.13, 'OB: %i | %s' % (fc_params['obId'],lab), relative=True, 
                 horizontalalignment='left')
   ax1.add_label(0.0,1.08, r'$\lambda_{fc}$: '+bk_lam, relative=True, 
                 horizontalalignment='left')
   
   # Display the observing date              
   #ax1.add_label(1.0,1.02, r'@ '+datetime.strftime(fcm_m.obsdate, '%d-%m-%Y %H:%M UTC%z'), relative=True, 
   #              horizontalalignment='right', fontsize=12)
   ax1.add_label(1.0,1.02, r'Obs. date: '+datetime.strftime(fcm_m.obsdate, '%Y-%m-%d'), relative=True, 
                 horizontalalignment='right', fontsize=12)
   
   #Finally include the version of fcmaker in there
   ax1.add_label(1.02,0.02, r'Created with fcmaker v%s'%(fcm_m.__version__), 
                 relative=True, 
                 horizontalalignment='left',verticalalignment='bottom',
                 fontsize=10, rotation=90)          

   # Save it all, both jpg for upload to p2, and pdf for nice quality.
   fn_out = os.path.join(fcm_m.plot_loc,fc_params['ob_name'].replace(' ','_')+'_'+ 
                         survey.replace(' ','-'))
                                       
   fig1.savefig(fn_out+'.jpg')
   if do_pdf:
      fig1.savefig(fn_out+'.pdf')
   if do_png:
      fig1.savefig(fn_out+'.png')
      
   plt.close()

   return fn_out+'.jpg'
