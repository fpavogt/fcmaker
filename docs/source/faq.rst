.. _faq:

FAQ 
====


**Q: If fcmaker is not an official ESO tool, will the fcmaker charts still be accepted in Phase 2?**

   **A**: fcmaker is developed by an ESO Fellow with duties on UT4 at the VLT, but it is *not* an official ESO tool. Still, the fcmaker charts match all the official ESO requirements for finding charts, and were fine-tuned to meet the specific needs of the telescope operators and night astronomers during operations. The hope is that the inclusion of fcmaker finding charts in OBs will thus remain tolerated. Even for those instruments for which the use of the GuideCam Tool may be compulsory, I would still encourage the use of
   fcmaker once the OB has been finalized.
   
   
**Q: How does fcmaker deal with proper motions?**
   
   **A**: fcmaker relies on ``astropy.coordinates.SkyCoord.appl_space_motion()`` routine to
   propagate proper motion between epochs. This routine assumes that the target is moving
   in a straight line (in space). It also requires the distance along the-line-of sight to
   the target. In comparison, the VLT propagates proper motions by assuming that the 
   target moves along a Great Circle on the sky. The error associated with this mismatch 
   will be negligible in most cases, particularly if user-provided coordinates are at 
   recent epochs.

   