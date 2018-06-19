.. _faq:

FAQ 
====   
   
**Q: How does fcmaker deal with proper motions?**
   
   **A**: fcmaker relies on ``astropy.coordinates.SkyCoord.appl_space_motion()`` routine to
   propagate proper motion between epochs. This routine assumes that the target is moving
   in a straight line (in space). It also requires the distance along the-line-of sight to
   the target. In comparison, the VLT propagates proper motions by assuming that the 
   target moves along a Great Circle on the sky. The error associated with this mismatch 
   will be negligible in most cases, particularly if user-provided coordinates are at 
   recent epochs.

   