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

#from . import fcmaker_muse as fcm_muse

'''
 fcmaker: a Python module to automatically create finding charts for ESO OBs in p2.\n
 Copyright (C) 2017,  F.P.A. Vogt
 --- oOo ---
 This file contains general functions for the fcmaker routines. 
 Created October 2017, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''

