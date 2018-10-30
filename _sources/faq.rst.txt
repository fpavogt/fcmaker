.. _faq:

FAQ 
====   

**Q**: What's with the versioning scheme ?

   **A**: Each version of fcmaker is flagged as ``vXXX.Y.Z``, starting from ``v103.0.0``. 
   ``XXX`` designs the observing period (in ESO convention) that fcmaker has been tuned to. 
   ``Y`` and ``Z`` indicate major and minor subversions, respectively. As a rule of thumb, 
   important updates that can justify to have users upgrade their version of the code will 
   increase ``Y``. In turn, minor upgrades (e.g. cosmetic changes with no impact on the 
   scientific content of the chart) will lead to an increase in ``Z``. 

**Q: How does fcmaker deal with proper motions?**
   
   **A**: fcmaker relies on ``astropy.coordinates.SkyCoord.appl_space_motion()`` routine to
   propagate proper motion between epochs. This routine assumes that the target is moving
   in a straight line (in space). It also requires the distance along the-line-of sight to
   the target. In comparison, the VLT propagates proper motions by assuming that the 
   target moves along a Great Circle on the sky. The error associated with this mismatch 
   will be negligible in most cases, particularly if user-provided coordinates are at 
   recent epochs.

   