# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------------------

import os
from dateutil import parser as dup
import astropy.units as u

'''
fcmaker: a Python module to automatically create finding charts for ESO OBs in p2.\n
Copyright (C) 2017-2018,  F.P.A. Vogt
--- oOo ---
This file contains some generic metadata used throughout the fcmaker module, including 
the version number, etc ...
Created October 2017, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''

# Define the version of fcmaker
__version__ = '0.3.0'

# Where are we located ?
fcm_dir = os.path.dirname(__file__)

# Specify the plotting style
fcm_usetex = True
fcm_plotstyle = os.path.join(fcm_dir,'mpl_styles','fcmaker_plots.mplstyle')

# Naming convention for downloading and storing stuff
plot_loc = os.path.join('fcm_plots') # Where to store the finding charts
data_loc = os.path.join('fcm_data') # Where to put the background images

# Radius of the downloaded image = max radius of VLT
bk_radius = 25*u.arcmin

# Set the pixel scale for the downloaded images
bk_pix = 0.75*u.arcsec

# Clear the SkyView cache ? Set to False, to avoid re-downloading images all the time.
clear_SkyView_cache = False

# Force the finding chart to be oriented N-S ? Requires Montage to be installed
set_North = True

# Default date (and time) at which the observations will take place.
obsdate = dup.parse("2018 07 01 00:00:00 UTC")

# VLT Guide star nominal magnitude range (in UCAC2 UCmag system)
gs_mag = [10.,14.]
# Maximum search radius for GS.
outer_GS_search = 11. * 60. # in arcsec

# minimum value of the GAIA propermotion stars to display
# set to -1 to never do this
min_abs_GAIA_pm = 100*u.mas/u.yr # in mas/yr

# The length of the tracks for high PM stars 
# Ideally, I would track it between the date of the finding chart and today ...
# But stupid SkyView does not return the exact DATE-OBS coming from the fits file ...
pm_track_time = 10 * u.yr

# For moving target, what time interval (+-) to show on the plot ?
ephem_range = 2*u.hour

# For moving targets, I need to assume a distance to propagate the proper motions between
# the closest ephemeris point and the actual observing time.
# Default should do a "reasonable job most of the times ..."
ephem_d = 1*u.au # Default 10 AU, because moving targets are likely to be in the Solar system
