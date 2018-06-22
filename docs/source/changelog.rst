.. _changelog:

.. |last-commit| image:: https://img.shields.io/github/last-commit/fpavogt/fcmaker.svg?colorB=e6c000
   :target: https://github.com/fpavogt/fcmaker

Changelog |last-commit|
=======================

.. todo:: 
   - (?) include an instrument-free mode
   - (?) add metadata to .jpg 
   - (!) formally validate the parallactic function
   - (?) add support for jitter in HAWKI (showing the max jitter area with a circle) 
   - (?) find a better way to display the allowed TT area for MUSE, e.g. shaded area with 
     ``shapley``
   - (?) make the obsdate an ``fc_params`` entry rather than a global variable
   - (?) for time critical OBs, get the time from the OB

v0.3.5 June 2018, F.P.A. Vogt
 - fixed bad bug when feeding no ``obid``

v0.3.4 June 2018, F.P.A. Vogt
 - fixed missing relsize package in mplstyle
 - added a check to make sure the user provides the wavelength of custom fits files
 - added section about custom FITS file to the doc.
 - fixed lack of `ins_mode` in new NFM OBs
 - added new keywords for command line option, including ``--p2uid``, ``--obsid``, 
   ``--plot-loc``, ``--data-loc``, ``--bk-image`` and ``--bk-lam``.

v0.3.3 June 2018, F.P.A. Vogt:
 - corrected the orientation of the MUSE WFM field (180 flip required)
 - added correct Cassegrain field-of-view of 7.4 arcmin in radius
 - formalized support for MUSE NFM, including ... 
 - the creation of mock bk_images from all the Gaia entries in the area.
 - fixed bug in radius of allowed GS displayed
 - fixed bug with chart tags (now displayed only when required)
 - used telescope coordinates from ESO website
 - added clear mark for slit PA in X-shooter charts

v0.3.2 May 2018, F.P.A. Vogt:
 - added \*.fits to .gitignore
 - added Pillow to the list of required packages (for direct jpg exports)
 - added a few more Exceptions to fool-proof stuff
 - bumped p2api version request to 0.92, to have the fix for ephemeris files
 - added do_parang keyword: by-default, hide the instrument field-of-view if a parallactic angle is required
 - updated Gaia DR2 article link in doc
 - added symbols for OBs with moving targets or parallactic angles
 - added \*.pdf and \*.jpg to .gitignore
 - added AO support for HAWKI
 - added RRM templates to MUSE and HAWKI

v0.3.1 May 2018, F.P.A. Vogt:
 - added XSHOOTER to the list of supported instruments
 - gave up on using the OBS-DATE keywords to draw the proper motion tracks for stars in
   the field. Always use fcm_m.pm_track_time instead.
    
v0.3.0 May 2018, F.P.A. Vogt:
 - added 2 functions to run fcmaker from within a script (make_fc and make fc_local)
 - restructured _main_.py and fcmaker.py as a result
 - replaced 'propagate_pm' with new 'SkyCoord.apply_space_motion()' function from Astropy 3.0
 - draw the proper motion vectors of the fastest stars in the field of view, using GAIA DR2.
 - for these, if OBS-DATE is in the fits header, then plot the pm line between obstime and 
   then. Else, plot as long as fcm_m.pm_track_time
 - started working on support for XSHOOTER
 - when DSS2 Red is not used for the zoomed-in view, still use it for the right-hand-side 
   plot
 - when reading a local file, only read the keywords that matter
 - made the use of Python-Latex and No-montage the default (safer for new users)
 - added support for moving targets with ephemeris files (on P2)

v0.2.1 January 2018, F.P.A. Vogt:
 - fixed a bad bug with the p2api import

v0.2.00 January 2018, F.P.A. Vogt:
 - remove local version of p2api in favor of pip one
 - initial Github+pypi release
 - fixed bug when no upload required (reply '0')
     
v0.1.48 November 2017, F.P.A. Vogt:
 - added Gallery page to docs, to show all local setup files, and some plots
 - fixed LaTeX bug when flagging bad telluric stars
 - added validity check of user Guide Stars (too close/far away ?)
 - added ability to export to png directly (helps with direct inclusion in docs)

v0.1.47 November 2017, F.P.A. Vogt:
 - added basic support for HAWKI, incl. Fast Phot acquisitions
 - added query to UCAC2 via Vizier, to show which Guide Stars are suitable
 - added check of TTS validity for MUSE AO (distance-wise), flagging the bad ones
 - added variable size of the right-hand-side plot, to show all offsets, even the very 
   large ones
 - add pypi badge to main page

v0.1.46 November 2017, F.P.A. Vogt:
 - fixed a bug when the length of offsets in smaller than ``noff``
 - fixed a bug when there is only one AO TTS defined.
 - added support for ``MUSE_wfm_cal_astrom`` and ``MUSE_wfm_cal_specphot``

v0.1.45 November 2017, F.P.A. Vogt:
 - updated doc with correct example chart

v0.1.43 November 2017, F.P.A. Vogt:
 - added ``Intended audience``, ``Topic`` and ``License`` flags to pypi release
 - implemented support for DETECTOR offsets in MUSE
 - added the ``--obsdate`` flag, to feed the date of the observation to fcmaker
 - added the ``--obsdate`` and ``__version__`` to the finding charts
 - added support of target proper motion of 1st order (assuming flat sky)
 - updated doc
   
v0.1.42 November 2017, F.P.A. Vogt:
 - in case of large blind offset, have a flexible zoom level in the left plot panel
 - added option to save to PDF in fcmaker_plots.make_fc() and __main__.py
 - fixed all docstrings in p2api.py
 - handle the orientation of custom background images by rotating the N-E arrows.
 - added a default "target" field for all MUSE WFM (AO) observations, because this is where
   the system first closes the AO loop (before applying any of the offset in the observing
   sequence). This allows the observer to check that the TTS are also valid in this position.
 - added legend to the charts
 - restructured the plotting to better separate instrument-specific elements from generic 
   ones. Created fcmaker_instrument_dispatch.py to that effect.
 - added support for OBs with multiple science templates
 - added the obId to chart

v0.1.41 November 2017, F.P.A. Vogt:
 - propagated the PA of the acquisition frame to the Science sequence (MUSE)
 - allowed to specify only 1 bk_image and bk_lams in automated mode

v0.1.40 October 2017, F.P.A. Vogt:
 - included doc in pypi package
 - updated doc

v0.1.26 October 2017, F.P.A. Vogt:
 - pre-release
 - initial doc assembled
 

 
  
 
