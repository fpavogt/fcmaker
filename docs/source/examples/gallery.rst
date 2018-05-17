.. _gallery:

Gallery
=======

All the fcmaker finding charts share common elements. They are:

 - **two plots per finding chart:** a zoom-in view on the left, and a global view on the right. 
 - **a right panel aimed at the telescope operator:** it always uses DSS2 Red images for the 
   background. 
   
   + Black circles, 11 arcmin in radius, mark the validity area for Guide Stars.
   + All suitable UCAC2 stars with 11<UCmag<15 within are tagged with red circles: these 
     stars are the ones that fit the nominal VLT GS specifications, and can be readily 
     selected by the Telescope Operator during the OB acquisition. Only the UCAC2 stars 
     compatible with all the OB offsets are tagged. 
   + If specified, the user-selected Guide Star is shown, and flagged accordingly if it 
     is unsuitable for any of the offset positions. 
   + Any star with large proper motion (where large is by default greater than 100 mas/year 
     according to GAIA DR2, set in ``fcm_m.min_abs_GAIA_pm``) is drawn, together with a 
     20 year look-back track (set in ``fcm_m.pm_track_time``) derived with the
     ``astropy.coordinates.SkyCoord.appl_space_motion()`` routine, assuming a default 
     distance of 100 pc (set in ``fcm_m.default_pm_d``).
     
 - **a left panel aimed at the night astronomer:** the background image is instrument  
   dependant, and can be chosen/provided by the user. 
   
   + The acquisition field is shown in bold purple. 
   + In case of moving targets (with an ephemeris file), the target position within
     2 hours (default, set in ``fcm_m.ephem_range``) from the ``obsdate`` parameter are 
     shown with red stars. Each time entry within the ephemeris file is shown individually.
     The target position at the chosen ``obsdate`` is propagated from the closest ephemeris 
     point in time using the ``astropy.coordinates.SkyCoord.appl_space_motion()`` routine,
     assuming a default distance of 1 AU (set in ``fcm_m.ephem_d``). The error associated 
     with this choice will remain negligible if there is an ephemeris point close-enough. 
   + In case of moving targets (with proper motions in arcsec/year), the target coordinates
     at ``obsdate`` are derived with the ``astropy.coordinates.SkyCoord.appl_space_motion()`` 
     routine, assuming a default distance of 100 pc (set in ``fcm_m.default_pm_d``). 
     **This is not how the VLT will compute the target coordinates, instead assuming that
     it moves along a Great Circle**. The implied error should however remain small for in
     most cases, and even more so if the epoch of the target is recent.

Below are some typical examples of fcmaker finding charts for all the supported instruments,
along with a brief description of each chart specificities.


MUSE
----

MUSE WFM-NOAO
.............

The MUSE finding charts in NOAO mode show the location of the *target* defined in the OB, 
the acquisition field, and the subsequent O and S fields. The minimum valid radius for 
telescope Guide Stars is 120 arcsec from any offset position. 

.. figure:: ./fcm_plots/MUSE_WFM_NOAO_pm_DSS2-Red.png
    :width: 750px
    :align: center
    :alt: MUSE NOAO

To recreate this example finding chart, download 
:download:`local_2_fcm.muse_wfm-noao <./local_2_fcm.muse_wfm-noao-pm>` and run::
   
   python -m fcmaker -l -f local_2_fcm.muse_wfm-noao-pm --do-png --systemtex


MUSE WFM-AO
...........

In addition to the elements of the MUSE NOAO charts, the AO charts display the validity 
area for the Tip-Tilt stars for all OB offsets (incl. the target itself, which is when the
AO loop is first closed). If one of the TTS falls outside of the suitable area for any of
the offset position, it is flagged with ``!``. 

.. figure:: ./fcm_plots/MUSE_WFM_AO_DSS2-Red.png
    :width: 750px
    :align: center
    :alt: MUSE NOAO

To recreate this example finding chart, download 
:download:`local_2_fcm.muse_wfm-ao <./local_2_fcm.muse_wfm-ao>` and run::
   
   python -m fcmaker -l -f local_2_fcm.muse_wfm-ao --do-png --systemtex


HAWKI
-----

HAWKI NOAO
..........

The HAWKI finding charts in NOAO mode show the acquisition field, and the subsequent O and 
S fields if the ``HAWKI_img_obs_GenericOffset`` template is used. The minimum valid radius 
for telescope Guide Stars is 240 arcsec from any offset position. 

.. figure:: ./fcm_plots/HAWKI_NOAO_2MASS-K.png
    :width: 750px
    :align: center
    :alt: HAWKI NOAO

To recreate this example finding chart, download 
:download:`local_2_fcm.hawki_noao <./local_2_fcm.hawki_noao>` and run::
   
   python -m fcmaker -l -f local_2_fcm.hawki_noao --do-png --systemtex


HAWKI NOAO FastPhot
...................

In addition to the elements of the HAWKI NOAO charts, the FastPhot charts display the  
detector windowed area. 

.. figure:: ./fcm_plots/HAWKI_NOAO_FPH_2MASS-K.png
    :width: 750px
    :align: center
    :alt: HAWKI NOAO FastPhot

To recreate this example finding chart, download 
:download:`local_2_fcm.hawki_fph <./local_2_fcm.hawki_fph>` and run::
   
   python -m fcmaker -l -f local_2_fcm.hawki_fph --do-png --systemtex

XSHOOTER
--------
The XSHOOTER finding charts show the field-of-view of the acquisition camera in bold purple.
Slit, IFU or acquisition camera field-of-views are then shown for each offset position, 
according to the selected observing templates. The minimum valid radius for telescope 
Guide Stars is 120 arcsec from any offset position.  If requested in the OB, fcmaker will
draw the XSHOOTER slit/IFU at the parallactic angle **at the time specified by** ``obsdate``.
In such cases when the finding charts becomes highly time-dependant, the observing time 
appears in red.

.. figure:: ./fcm_plots/XSHOOTER_OB_EPHEM_DSS2-Red.png
    :width: 750px
    :align: center
    :alt: XSHOOTER EPHEM

To recreate this example finding chart, download 
:download:`local_2_fcm.xshooter_ephem <./local_2_fcm.xshooter_ephem>` and run::
   
   python -m fcmaker -l -f local_2_fcm.xshooter_ephem --do-png --systemtex --obsdate 2018-05-15 08:23:00 UTC



