.. _examples:
.. highlight:: python

Usage instructions
==================

Case 1: manual creation of a finding chart via p2
-----------------------------------------------------

In its most basic mode, fcmaker connects to `p2 <http://www.eso.org/p2>`_ and creates the finding chart for a given 
OB Id. To do so, simply type in a  terminal::

   python -m fcmaker

or from within a Python shell::

   
   >>> run -m fcmaker
   
You will be prompted for your `p2 <http://www.eso.org/p2>`_ username, password, and the 
Id of the observing block to process. fcmaker will create two folders ``fcm_data`` and 
``fcm_plots`` at your current location, where it will store the background image and the 
finding charts. Once the finding charts have been generated, you can choose to attach the 
newly created finding chart to the OB on `p2 <http://www.eso.org/p2>`_, or not.

Case 2a: semi-automatic creation of finding charts via p2 (incl. upload)
------------------------------------------------------------------------

If you have a lot of finding charts to create, fcmaker allows for a semi-batch processing.
First, create a text file :download:`p2_2_fcm.params <./examples/p2_2_fcm.params>` (the actual filename is flexible) with the following structure:

.. literalinclude:: examples/p2_2_fcm.params

Then, feed it to fcmaker with the ``-f`` flag::

   python -m fcmaker -f p2_2_fcm.params

or from within a Python shell::
   
   >>> run -m fcmaker -- -f p2_2_fcm.params
   

**Note the extra** ``--`` **before listing any of the** `flags`_ **with the** ``run`` **command!**

In doing so, fcmaker will connect to `p2 <http://www.eso.org/p2>`_ and process all the 
OBs listed. Note that for each finding chart, you will still need to manually specify 
whether you want to upload it to `p2 <http://www.eso.org/p2>`_ (or not). 

Case 2b: fully automatic creation of finding charts from p2 (no upload)
-----------------------------------------------------------------------

You can fully automate the creation of many finding charts if you include the 
``--no-upload`` flag. In that case, fcmaker will never attempt to upload anything to 
`p2 <http://www.eso.org/p2>`_, and only save them locally::

   python -m fcmaker --no-upload -f p2_2_fcm.params

Case 3: targets with proper motions
-----------------------------------

In case of large proper motions of the target, one can provide the expected year, (month, 
day ...!) of the observation to create an accurate finding chart, using the ``--obsdate`` 
flag. For example::

   python -m fcmaker -f randomfilename --obsdate 2018
   python -m fcmaker -f randomfilename --obsdate 2018-05
   python -m fcmaker -f randomfilename --obsdate 2018-05-17
   
   
Case 4: creation of finding charts locally (without p2)
-------------------------------------------------------

fcmaker can also create finding charts *locally*, without the need to have an OB present 
on `p2 <http://www.eso.org/p2>`_ first. To do so, one first needs to specify the basic OB 
parameters in a text file ``whichevername`` (in essence, a stripped-down version of a 
full OBX file):
   
   - :download:`local_2_fcm.params.muse  <./examples/local_2_fcm.muse_wfm-ao>`
   - :download:`local_2_fcm.params.hawki <./examples/local_2_fcm.hawki_noao>`

The file is then fed to fcmaker with the ``-f`` flag, together with the ``--local`` flag to 
indicate that it is a *local run*::

   python -m fcmaker --local -f local_2_fcm.muse_wfm-ao
   python -m fcmaker --local -f local_2_fcm.hawki_noao
      
fcmaker will create the associated finding chart, store it where specified, and exit.

.. _flags:

The fcmaker flags
-----------------

A series of flags allow to fine-tune the way fcmaker works. They are:
 
 * ``--help,-h`` : prints the basic help
 * ``--version`` : prints the fcmaker version
 * ``--no-montage`` : will by-pass the use of Montage. Finding charts will rotated according to the input FITS file. [no encouraged]
 * ``--no-systemtex`` : will use the matplotib LaTeX, rather than the system-wide LaTeX installation. Finding charts will be less pretty.
 * ``--no-upload`` : tell fcmaker never to try to upload the finding chart(s) to p2.
 * ``--do-pdf`` : tell fcmaker to save a .pdf version of the finding chart (in addition to the jpg one).
 * ``--do-png`` : tell fcmaker to save a .png version of the finding chart (in addition to the jpg one).
 * ``--obsdate`` : allow to specify the observing date (and time), in case of large proper motions.
 * ``--clear-SkyView-cache`` : does as it says. Note that the SkyView cache is *not* the same as the ``fcm_data`` folder in which fcmaker stores the downloaded FITS files.


The background images
----------------------

fcmaker relies on ``astroquery.skyview`` to download background images for the finding 
charts (if no local FITS file is provided by the user). To display the full list of 
surveys available, type in a Python shell::
 
   from astroquery.skyview import SkyView
   SkyView.survey_dict['overlay_blue']

The default background survey images for the instruments supported by fcmaker are as follows:
   
   * **MUSE**: ``DSS2 Red``
   * **HAWKI**: ``2MASS-J``, ``2MASS-H`` and ``2MASS-K`` for filters ``J``, ``H``, ``K``, and ``2MASS-H`` for all other filters

.. note::
   
   From an operational perspective, when using Skyview images, I would strongly recommend 
   to only use ``DSS2 Red`` or ``SDSSr`` for optical instruments, and ``2MASS-J``, 
   ``2MASS-H`` or ``2MASS-K`` for IR instruments.

