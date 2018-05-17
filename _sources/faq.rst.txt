.. _faq:

FAQ 
====


**Q: If fcmaker is not an official ESO tool, will the fcmaker charts still be accepted in Phase 2?**

   **A**: fcmaker is developed by an ESO Fellow with duties on UT4 at the VLT, but it is *not* an official ESO tool. Still, the fcmaker charts match all the official ESO requirements for finding charts, and were fine-tuned to meet the specific needs of the telescope operators and night astronomers during operations. The hope is that the inclusion of fcmaker finding charts in OBs will thus remain tolerated. Even for those instruments for which the use of the GuideCam Tool may be compulsory, I would still encourage the use of
   fcmaker once the OB has been finalized.
   
   
**Q: How does fcmaker deal with proper motions?**
   
   **A**: fcmaker relies on ``astropy.coordinates.SkyCoord.appl_space_motion()`` routine to
   propagate proper motion between epochs. This routine assumes that the target is moving
   in a straight line. It also requires a distance along-the-line-of sight. In comparison, 
   the VLT propagate proper motions by assuming that the target moves along a Great Circle
   on the sky, and thus does **not** require a line-of-sight distance. 
   
   By default, fcmaker assumes a distance of 1 AU (``fcm_m.ephem_d``) for targets with an 
   ephemeris files, which by-and-large are found within the Solar System. For targets with
   proper motions defined as arcsec/year (either GAIA DR2 entries or the target set by the 
   user), fcmaker assumes a default distance of 100 pc (``fcm_m.default_pm_d``). The error
   associated with these choices will remain negligible in most cases, particularly if 
   user-provided coordinates are at recent epochs.

   