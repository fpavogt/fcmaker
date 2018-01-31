.. _changelog:

Changelog
==========

.. todo:: 
   - add support for moving target (i.e. deal with ephemeris files)
   - add support for jitter in HAWKI (showing the max jitter area with a circle)
   - deal with proper motion in a more accurate way than in `propagate_pm`
     (once ``astropy v3.0`` is out ?)
   - draw the proper motion vectors of the fastest stars in the field of view
     (once GAIA v2 will is out ?)   
   - find a better way to display the allowed TT area for MUSE ? 
     e.g. shaded area ?

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
 

 
  
 
