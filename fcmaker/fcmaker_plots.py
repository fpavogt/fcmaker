# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------------------

import os
import numpy as np

import warnings # for aplpy displaying nasty packages
from datetime import datetime
from dateutil import parser as dup
import pytz
import copy

# Import plotting packages
from matplotlib import pylab as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines

with warnings.catch_warnings():
   warnings.simplefilter("ignore")
   import aplpy

# Some astro stuff
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy import units as u
import astropy.io.fits as fits
from astropy.constants import c
from astropy.time import Time
from astroquery.skyview import SkyView
from astroquery.vizier import Vizier
from astroquery.gaia import Gaia

from astropy.wcs import WCS

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
                  #head_length = 562 * arrow_width.to(u.deg).value,
                  #head_width = 422 * arrow_width.to(u.deg).value, 
                  color='k')
   
   #print(arrow_width.to(u.deg).value)
   
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
def get_bk_image(bk_image, bk_lam, center, radius, fc_params):
   ''' 
   A function that fetches the background image from the web, or locally if specified
   by the user.
   
   Args:
      bk_image: (string) defines whether to use a local image (FITS filename) or a SkyView 
                image (survey name). 
      bk_lam: (string)
      center: a SkyCoord entry defining the center of the image
      radius: the radius of the image ... for now, only used when creating a pseudo-Gaia image
      fc_params: the parameters of the finding chart (a dictionary)
   Returns:
      The absolute path of the background image (possibly after download), and the survey name.
   
   .. note::
     
      For local images, the code expect them to be stored in `fcm_m.data_loc` folder.
      
   '''
   
   if not(bk_image in [None,'None']):
      # Has the user provided his/her own FITS file ?  
      if not(bk_image in SkyView.survey_dict['overlay_blue']) and not(bk_image=='Gaia'):
         # Very well, I guess I have a custom FITS file.
         fn_bk_image = os.path.join(fcm_m.data_loc, bk_image)
         
         #survey = bk_image.split('.')[0] # assume a filename.fits ...
         survey = 'from user'
         # Here, make sure the user also provides the wavelength of the image
         if bk_lam in [None,'None']:
            raise Exception('Ouch! Please specify the wavelength of the background image '+ 
                             'with the bk_lam(s) parameter.')
         
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
   # If not, get one
   if fcm_m.clear_SkyView_cache or not(os.path.isfile(fn_bk_image)):
   
      # Specify FITS filename
      fn_bk_image = os.path.join(fcm_m.data_loc, fc_params['ob_name'].replace(' ','_')+'_' + 
                                 survey.replace(' ','-') + 
                                 '.fits')
      
      
      if survey == 'Gaia': # Create a home-made Gaia FITS image
         
         (sampling,fwhm) = fcm_id.get_gaia_image_params(fc_params)
         
         make_gaia_image(center, 1.1*2*radius*u.arcsec, # 1.1 because plot is rectangular, 2x to get diameter
                         sampling, fwhm, 
                         obsdate=fcm_m.obsdate, fn = fn_bk_image)
         
      else: # Use SkyView to get the image 
         print('   Downloading the %s background image with SkyView ...' % (survey))
         path = SkyView.get_images(position=center, survey=[survey], 
                                   radius=fcm_m.bk_radius,
                                   pixels='%i' %((fcm_m.bk_radius/fcm_m.bk_pix).to(u.arcsec/u.arcsec).value),
                                   )
         # Use a direct download link for DSS2
         # use urllib and urlretrieve
         # https://archive.stsci.edu/cgi-bin/dss_search?v=poss2ukstu_ir&r=0+0+0+&d=-0+0+1&e=J2000&h=15.0&w=15.0&f=gif&c=none&fov=NONE&v3=
         # https://archive.stsci.edu/cgi-bin/dss_search?v=poss2ukstu_ir&r=0+0+0+&d=0+0+0&e=J2000&h=15.0&w=15.0&f=fits&c=none&fov=NONE&v3=                          
                              
         if len(path)==0:
            raise Exception('No image downloaded via SkyView.get_images!')
      
         # Save it locally for later
         outfits = fits.HDUList([path[0][0]])
         outfits.writeto(fn_bk_image, output_verify='fix', overwrite=True)

   if bk_lam in [None,'None']:
      bk_lam = get_bk_image_lam(fn_bk_image, fc_params)

     
   return (fn_bk_image, survey, bk_lam)
   
# ----------------------------------------------------------------------------------------
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
   # Fixes issues #4: https://github.com/fpavogt/fcmaker/issues/4
   freqs = [ l.split(' THz')[0].split(' ')[-1].split('-') for l in header['COMMENT'] if ('Bandpass' in l) and ('THz' in l)][0]
   # WARNING HERE: I assume the bandpass is in THz in all cases !!!!!
   if len(freqs) !=2:
      raise Exception('Ouch! Something went wrong when extracting the frequency of the background image!')
      
   freqs = [np.int(c.value/(np.float(freq)*1.e12)*1.e9) for freq in freqs]
   bk_lam = r'%i - %i nm' % (freqs[1],freqs[0])

   return bk_lam

# ----------------------------------------------------------------------------------------
def make_gaia_image(skycoord, fov, sampling, fwhm,
                    obsdate = dup.parse("2018 07 01 00:00:00 UTC"), 
                    fn = os.path.join('.','pseudo_gaia.fits')):
   ''' 
   Given a coordinate, use the Gaia catalogue to create a mock image of the sky.
   
   Args:
      skycoord: the astropy SkyCoord of the center of the image
      fov: size of the image to create, in astropy.units
      sampling: pixel size of the image to create, in astropy.units
      fwhm: fwhm of the fake stars in the image, in astropy.units
      obsdate: date of the observation (to propagate the motion of the stars) as datetime.datetime object 
      fn: relative path + FITS filename
   Returns:
      The FITS filename
   '''

   # Ok, assemble an array of the suitable size
   nx = int(np.floor((fov/sampling).value))
   im = np.zeros((nx,nx))

   # Query GAIA
   print('   Querying GAIA to create a pseudo-image of the sky ...')
   # Make it a sync search, because I don't think I need more than 2000 targets ...
   
   j = Gaia.cone_search(skycoord, fov, verbose=False)
   r = j.get_results()
   nstars = len(r)
   
   # TODO: well, what if I do need more than 2000 targets ?
   # Issue a warning for now ...
   if nstars == 2000:
      warnings.warn(' Reached query limit of 2k stars. Gaia image will not be complete.')

   # Create a very basic header for the FITS file
   hdr = fits.Header()
   hdr.set('CRPIX1', nx/2, 'X reference pixel')
   hdr.set('CRPIX2', nx/2, 'Y reference pixel')
   hdr.set('CRVAL1', skycoord.ra.deg, 'Reference longitude')
   hdr.set('CRVAL2', skycoord.dec.deg, 'Reference latitude')
   hdr.set('CTYPE1', 'RA---TAN', 'Coordinates -- projection')
   hdr.set('CTYPE2', 'DEC--TAN', 'Coordinates -- projection')
   hdr.set('CDELT1', -sampling.to(u.deg).value, 'X scale')# mind the sign !
   hdr.set('CDELT2', sampling.to(u.deg).value, 'Y scale')
   hdr.set('RADESYS', 'ICRS', 'Coordinate system') 
   hdr.set('EQUINOX', 2000.0, 'Epoch of the equinox')
   hdr['COMMENT'] = 'Artificial sky image reconstructed from the Gaia catalogue.'
   hdr['COMMENT'] = 'Number of stars: %i' % (nstars)
   hdr['COMMENT'] = 'FWHM: %i mas' % (int(fwhm.to(u.mas).value))
   hdr['COMMENT'] = 'Bandpass:    285.5-908 THz '
   hdr['COMMENT'] = 'Epoch: %s' % (obsdate.strftime('%Y-%h-%d %H:%M:%S UTC'))
   hdr['COMMENT'] = 'Created %s with fcmaker v%s' % (datetime.utcnow().strftime('%Y-%h-%d %H:%M:%S UTC'),
                                                        fcm_m.__version__)
   hdr['COMMENT'] = 'See http://fpavogt.github.io/fcmaker for details.'
   
   # Save the empty FITS file for now
   hdu = fits.PrimaryHDU(im, header=hdr)
   hdul = fits.HDUList([hdu])
   hdul.writeto(fn, overwrite=True)
   
   # Re-extract the WCS info
   w = WCS(fn)

   # Prepare some variables  
   x, y = np.meshgrid(np.arange(0,nx,1), np.arange(0,nx,1))
   sigma = (fwhm/sampling).value/(2*np.sqrt(2*np.log(2)))

   # For each star, add a Gaussian of the suitable size
   for s in range(nstars):
      # First, create a SkyCoord entry
      star = SkyCoord(ra=r['ra'][s], dec=r['dec'][s], 
                 obstime=Time(r['ref_epoch'][s], format='decimalyear'),
                 frame = 'icrs', unit=(u.deg, u.deg), 
                 pm_ra_cosdec = r['pmra'][s]*u.mas/u.yr,
                 pm_dec = r['pmdec'][s]*u.mas/u.yr,
                 # I must specify a generic distance to the target,
                 # if I want to later on propagate the proper motions
                 distance = fcm_m.default_pm_d,  
                         )
   
      # Then, propagate the star at the time of the observation
      star_now = star.apply_space_motion(new_obstime = Time(obsdate))
      
      # Get the image coordinate of the star
      star_coords = w.wcs_world2pix([star_now.ra.deg],[star_now.dec.deg],0)
      
      # The distance of all pixels to the star
      d = np.sqrt((x-star_coords[0])**2 + (y-star_coords[1])**2)
      
      # Add the star to the array ... scale it as a function of its flux
      im += r['phot_g_mean_flux'][s]*np.exp(-( (d)**2 / ( 2.0 * (fwhm/sampling).value**2 ) ) )
   
   
   # Save the final FITS file   
   hdu = fits.PrimaryHDU(im, header=hdr)
   hdul = fits.HDUList([hdu])
   hdul.writeto(fn, overwrite=True)   
   
   return fn

# ----------------------------------------------------------------------------------------
def draw_fc(fc_params, bk_image = None, bk_lam = None, do_pdf = False, do_png = False):
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
   
   # Get the radius and center of the charts
   (left_radius, right_radius) = fcm_id.get_chart_radius(fc_params)
   (left_center, right_center) = fcm_id.get_chart_center(fc_params)
   
   # Get the background image in place
   (fn_bk_image_L, survey_L, bk_lam_L) = get_bk_image(bk_image, bk_lam, left_center, left_radius, fc_params)
   
   # If this is not the DSS2 Red, then download this as well for the right-hand-side plot
   if not(survey_L == 'DSS2 Red'):
      (fn_bk_image_R, survey_R, bk_lam_R) = get_bk_image('DSS2 Red', None, right_center, right_radius, fc_params)
   else: 
      (fn_bk_image_R, survey_R, bk_lam_R) = copy.deepcopy((fn_bk_image_L, survey_L, bk_lam_L))
   
   # Start the plotting
   plt.close(1)
   fig1 = plt.figure(1, figsize=(14.17,7)) #14.17 inches, at 50% = 1 full page plot in A&A. 

   ax1 = aplpy.FITSFigure(fn_bk_image_L, figure=fig1, north=fcm_m.set_North,
                          subplot=[0.12,0.12,0.35,0.76])
   ax1.show_grayscale(invert = True, stretch='linear', pmin = fcm_id.get_pmin(survey_L))

   ax2 = aplpy.FITSFigure(fn_bk_image_R, figure=fig1, north=fcm_m.set_North, 
                          subplot=[0.59,0.12,0.38,0.8])
   ax2.show_grayscale(invert = True, stretch='linear', pmin = fcm_id.get_pmin(survey_R))
   
   ax1.recenter(left_center.ra,left_center.dec,radius=left_radius/3600.)
   ax2.recenter(right_center.ra,right_center.dec,radius=right_radius/3600.) 
   
   # I will be adding stuff to the left-hand-side plot ... start collecting them
   ax1_legend_items = []
   
   # Do I have the epoch of observation for the background image ?
   '''
   try:
      bk_obsdate = fits.getval(fn_bk_image, 'DATE-OBS')
      
      # If not specified, assume UTC time zone for bk image   
      if bk_obsdate.tzinfo is None:
      #   #warnings.warn(' "--obsdate" timezone not specified. Assuming UTC.')
         bk_obsdate = bk_obsdate.replace(tzinfo=pytz.utc)
      
      nowm_time = Time(bk_obsdate)
      pm_track_time = (fcm_obsdate - bk_obsdate).total_seconds()*u.s
      
   except:
      # If not, just shows the default length set in fcmaker_metadata
      nowm_time = Time(fcm_m.obsdate) - fcm_m.pm_track_time
      pm_track_time = fcm_m.pm_track_time
   '''
   
   # Just keep things simple. Same lookback time for all charts.
   nowm_time = Time(fcm_m.obsdate) - fcm_m.pm_track_time
   #pm_track_time = fcm_m.pm_track_time
   
   # If yes, then let's show where are the fastest stars moved from/to.
   #if do_GAIA_pm and (fcm_m.min_abs_GAIA_pm >=0):
   if (fcm_m.min_abs_GAIA_pm >=0):
         
      # Query GAIA
      print('   Querying GAIA DR2 to look for high proper motion stars in the area ...')
   
      # Make it a sync search, because I don't think I need more than 2000 targets ...
      j = Gaia.cone_search(left_center, right_radius*u.arcsec, verbose=False)
      r = j.get_results()
   
      selected = np.sqrt(r['pmra']**2+r['pmdec']**2)*u.mas/u.yr > fcm_m.min_abs_GAIA_pm
   
      # Show the fastest stars
      past_tracks = []
      future_tracks = []
      
      # Need to propagate their coords to the fc epoch and the obstime
      for s in range(len(r['ra'][selected])):
      
         star = SkyCoord(ra=r['ra'][selected][s], dec=r['dec'][selected][s], 
                         obstime=Time(r['ref_epoch'][selected][s], format='decimalyear'),
                         frame = 'icrs', unit=(u.deg, u.deg), 
                         pm_ra_cosdec = r['pmra'][selected][s]*u.mas/u.yr,
                         pm_dec = r['pmdec'][selected][s]*u.mas/u.yr,
                         # I must specify a generic distance to the target,
                         # if I want to later on propagate the proper motions
                         distance=fcm_m.default_pm_d,  
                         )
      
         now = star.apply_space_motion(new_obstime = Time(fcm_m.obsdate))    
         nowm = star.apply_space_motion(new_obstime = nowm_time)
      
         past_tracks += [np.array([[nowm.ra.deg,now.ra.deg],[nowm.dec.deg,now.dec.deg]])]
         
         for ax in [ax1,ax2]:
            ax.show_markers([now.ra.deg],[now.dec.deg],marker='.',color='crimson', 
                             facecolor='crimson', edgecolor='crimson')
         
      #if len(r['ra'][selected])>0:
      #   # Prepare a dedicated legend entry
      #   ax1_legend_items += [mlines.Line2D([], [],color='crimson',
      #                         markerfacecolor='crimson',
      #                         markeredgecolor='crimson', 
      #                         linestyle='-',
      #                         linewidth=0.75,
      #                         marker='.',
      #                         #markersize=10, 
      #                         label='PM* (track$=-$%.1f yr)' % (pm_track_time.to(u.yr).value)) ]
      for ax in [ax1,ax2]:
         ax.show_lines(past_tracks,color='crimson',linewidth=0.75, linestyle = '-')
   
   # Query UCAC2 via Vizier over the finding chart area, to look for suitable Guide Stars
   print('   Querying UCAC2 via Vizier to look for possible Guide Stars ...')
   Vizier.ROW_LIMIT = 10000
   gs_table = Vizier.query_region(right_center, radius=right_radius*u.arcsec,
                                  inner_radius = fcm_id.get_inner_GS_search(fc_params)*u.arcsec,
                                  catalog="UCAC2")
   gs_table = gs_table[gs_table.keys()[0]]
   
   # Turn the table into a list of SkyCoord, that I can clean as I go.
   # Only select guide stars in the nominal GS mag range
   gs_list = [SkyCoord( ra=line['RAJ2000'], dec=line['DEJ2000'], obstime = Time('J2000'), 
                        equinox = 'J2000', frame = 'icrs', unit=(u.deg, u.deg), 
                        pm_ra_cosdec = line['pmRA']/1000.*u.mas/u.yr, pm_dec = line['pmDE']/1000.*u.mas/u.yr,
                        # Assume a fixed distance, so that I can then propagate proper motions
                        distance = 100.*u.pc).apply_space_motion(new_obstime = Time(fcm_m.obsdate))
              for line in gs_table if fcm_m.gs_mag[0] <= line['UCmag'] <= fcm_m.gs_mag[1]]
   
   
   # Here, I will show all the ephemeris points I have prepared, ahead of the 
   # observations (If I have any)
   
   if len(fc_params['ephem_points_past']) > 0:
      ax1.show_markers([item.ra.deg for item in fc_params['ephem_points_past']],
                       [item.dec.deg for item in fc_params['ephem_points_past']], 
                       marker='*', color='crimson', s=100,
                       facecolor='none', edgecolor='crimson')
      
      # Prepare a dedicated legend entry
      ax1_legend_items += [mlines.Line2D([], [],color='crimson',
                               markerfacecolor='none',
                               markeredgecolor='crimson', 
                               linestyle='',
                               #linewidth=0.75,
                               marker='*',
                               #markersize=10, 
                               label='Target ($\Delta T=-%.0f$ min)' % 
                                     (fc_params['ephem_past_delta'].total_seconds()/60.))]
                                     
   
   if len(fc_params['ephem_points_future']) > 0: 
      ax1.show_markers([item.ra.deg for item in fc_params['ephem_points_future']],
                       [item.dec.deg for item in fc_params['ephem_points_future']], 
                       marker='*', color='crimson', s=100,
                       facecolor='crimson', edgecolor='crimson')                 
      
      # Prepare a dedicated legend entry
      ax1_legend_items += [mlines.Line2D([], [],color='crimson',
                               markerfacecolor='crimson',
                               markeredgecolor='crimson', 
                               linestyle='',
                               #linewidth=0.75,
                               marker='*',
                               #markersize=10, 
                               label='Target ($\Delta T=+%.0f$ min)' % 
                                     (fc_params['ephem_future_delta'].total_seconds()/60.))]
             
   # Show the observation footprint
   for f in fields:
      
      # Call the function that will plot all the important stuff for this field.
      fcm_id.plot_field(ax1, ax2, fc_params, fields[f])
      
      # Keep all the possible guide stars in the area.
      gs_list = [star for star in gs_list if 
                 (fields[f][2].separation(star) > (fcm_id.get_inner_GS_search(fc_params)/3600.*u.deg)) 
                 and (fields[f][2].separation(star) < (fcm_id.get_GS_outer_radius(fc_params)/3600.*u.deg))]    
                   

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
   add_orient(ax2, right_center, radius = (right_radius*0.82)*u.arcsec, 
              arrow_width = (right_radius*0.82)/4.5 * u.arcsec, usetex = fcm_m.fcm_usetex)
   add_orient(ax1, left_center, radius = (left_radius*0.82)*u.arcsec, 
              arrow_width = (left_radius*0.82)/4.5 * u.arcsec, usetex = fcm_m.fcm_usetex)
   
    
   # Add a scale bar
   (scl,scll) = fcm_id.get_scalebar(fc_params['inst'], ins_mode = fc_params['ins_mode'])
   
   ax1.add_scalebar(scl) # Length in degrees
   ax1.scalebar.show(scl, label=scll, corner = 'bottom left', 
                     color = 'k', frame=1, fontsize=12)
   
   scl2 = np.floor(right_radius/60/6)/60.
   scll2 = r'%.0f$^{\prime}$' % (scl2*60)
   ax2.add_scalebar(scl2)
   ax2.scalebar.show(scl2, label=scll2, corner = 'bottom left', 
                     color = 'k', frame=1, fontsize=12)
   
   
   for ax in [ax1,ax2]:
      ax.scalebar.set_linewidth(2)
   
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
   ax1.add_label(0.0,1.11, 'Run ID: '+fc_params['prog_id']+' | '+fc_params['pi'], 
                 relative=True, horizontalalignment='left', size=14)
   
   # Fix some bugs for anyone not using LaTeX ... sigh ...
   if fcm_m.fcm_usetex:
      lab = fc_params['ob_name'].replace('_','\_')
   else:
      lab = fc_params['ob_name']             
    
          
   ax1.add_label(0.0,1.06, 'OB: %i | %s' % (fc_params['ob_id'],lab), relative=True, 
                 horizontalalignment='left', size=14)
   #ax1.add_label(0.0,1.08, r'$\lambda_{fc}$: %s (%s)' % (bk_lam_L, survey_L), relative=True, 
   #              horizontalalignment='left')
   
   # Display the Wavelength of the plots     
   ax1.add_label(0.02,0.965, r'%s (%s)' % (bk_lam_L, survey_L), relative=True, fontsize=12,
                 horizontalalignment='left',verticalalignment='center', bbox=dict(edgecolor='w', facecolor='w', alpha=0.85))
   ax2.add_label(0.02,0.965, r'%s (%s)' % (bk_lam_R, survey_R), relative=True, fontsize=12,
                 horizontalalignment='left',verticalalignment='center', bbox=dict(edgecolor='w', facecolor='w', alpha=0.85))
   
   
   # Add a legend (if warranted) for the left plot
   if len(ax1_legend_items) >0:
      ax1._ax1.legend(handles=ax1_legend_items,
                     bbox_to_anchor=(-0.03, 0.82, 0.4, .1), #loc='lower right',
                     ncol=1, #mode="expand", 
                     borderaxespad=0., fontsize=10, borderpad=0.3, 
                     handletextpad=0., handlelength=2.0)
   
   # Start keeping track of any tags I need to show
   tag_string = r' '
   
   # Show the obsdate
   ax1.add_label(1.0,1.02, r'Obs. date: '+datetime.strftime(fcm_m.obsdate, '%Y-%m-%d %H:%M %Z'), 
                 relative=True, color='k',
                 horizontalalignment='right', fontsize=11) 

   # Show the OB tags
   if 'moving_target' in fc_params['tags']:
      tag_string += '$\leadsto$ '
      
   if 'parallactic_angle' in fc_params['tags']:
      tag_string += '$\measuredangle$ '
   
   if len(tag_string)>1 : # only show the tag if it has a non-zero length
      ax1.add_label(-0.18,1.08, tag_string, 
                    relative=True, color='k',
                    horizontalalignment='center', verticalalignment='center', fontsize=30, 
                    bbox=dict(edgecolor='k', facecolor='lightsalmon', alpha=1, linewidth=1, 
                              boxstyle="sawtooth,pad=0.2,tooth_size=0.075"),
                    ) 

   # Finally include the version of fcmaker in there
   ax1.add_label(1.01,0.00, r'Created with fcmaker v%s'%(fcm_m.__version__), 
                 relative=True, 
                 horizontalalignment='left',verticalalignment='bottom',
                 fontsize=10, rotation=90)          

   # Save it all, both jpg for upload to p2, and pdf for nice quality.
   fn_out = os.path.join(fcm_m.plot_loc,fc_params['ob_name'].replace(' ','_')+'_'+ 
                         survey_L.replace(' ','-'))
                                       
   fig1.savefig(fn_out+'.jpg')
   if do_pdf:
      fig1.savefig(fn_out+'.pdf')
   if do_png:
      fig1.savefig(fn_out+'.png')
      
   plt.close()

   return fn_out+'.jpg'
