.. _examples:
.. highlight:: python

Usage instructions
==================

There are two main ways to use fcmaker. You can either execute the module as a script, i.e.::

   python -m fcmaker
   
or you can import the module and execute it within a script, i.e.::

   >>> import fcmaker
   >>> fcmaker.make_fc( ... ) # or fcmaker.make_fc_local ( ... )
   
When running fcmaker as a script, any arguments you feed it is pretty much sent to the 
underlying functions ``make_fc()`` and ``make_fc_local()``. For simplicity, this page
only discusses how to run the entire module as a script, which ought to be slightly more 
friendly to users not (yet!) familiar with Python. Still, the hope is that after reading 
this page, the use of the functions ``make_fc()`` and ``make_fc_local()`` shouldn't be
too mysterious. Also, don't forget the built-in help::

   >>> import fcmaker
   >>> help(fcmaker.make_fc)
   

Case 1: manual creation of a finding chart via p2
-----------------------------------------------------

In its most basic mode, fcmaker connects to `p2 <http://www.eso.org/p2>`_ and creates the finding chart for a given 
OB Id. To do so, simply type in a  terminal::

   python -m fcmaker

or from within a Python shell::
 
   >>> run -m fcmaker
   
You will be prompted for your `p2 <http://www.eso.org/p2>`_ username, password, and the 
ID of the observing block to process. fcmaker will create two folders ``fcm_data`` and 
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
day, ...!) of the observation to create an accurate finding chart, using the ``--obsdate`` 
flag. This also works for OBs that have Ephemeris files. For example::

   python -m fcmaker -f randomfilename --obsdate 2018
   python -m fcmaker -f randomfilename --obsdate 2018-05
   python -m fcmaker -f randomfilename --obsdate 2018-05-17 14:34:57 UTC
   
   
Case 4: creation of finding charts locally (without p2)
-------------------------------------------------------

fcmaker can also create finding charts *locally*, without the need to have an OB present 
on `p2 <http://www.eso.org/p2>`_ first (fcmaker will still require an internet connection, 
though!). To do so, one first needs to specify the basic OB parameters in a text file  
(in essence, a stripped-down version of a full OBX file). Here are
templates for all supported instruments:
   
   - :download:`local_2_fcm.params.muse  <./examples/local_2_fcm.muse>`
   - :download:`local_2_fcm.params.hawki <./examples/local_2_fcm.hawki>`
   - :download:`local_2_fcm.params.xshooter <./examples/local_2_fcm.xshooter>`
   
The file is then fed to fcmaker with the ``-f`` flag, together with the ``--local`` flag to 
indicate that it is a *local run*::

   python -m fcmaker --local -f local_2_fcm.muse
   python -m fcmaker --local -f local_2_fcm.hawki
   python -m fcmaker --local -f local_2_fcm.xshooter
   
fcmaker will create the associated finding chart, store it where specified, and exit.

.. _flags:

The fcmaker flags
-----------------

A series of flags allow to fine-tune the way fcmaker works. They are:
 
 * ``--help,-h`` : prints the basic help
 * ``--version`` : prints the fcmaker version
 * ``--montage`` : will use Montage, to force the finding charts to be rotated North (soon to be retired).
 * ``--systemtex`` : will use the system-wide LaTeX, rather than the matplotlib one. Finding charts will be prettier. [encouraged!]
 * ``--no-upload`` : tell fcmaker never to try to upload the finding chart(s) to p2.
 * ``--do-pdf`` : tell fcmaker to save a .pdf version of the finding chart (in addition to the jpg one).
 * ``--do-png`` : tell fcmaker to save a .png version of the finding chart (in addition to the jpg one).
 * ``--obsdate`` : allows to specify the observing date (and time), in case of large proper motions or moving target.
 * ``--clear-SkyView-cache`` : does as it says. Note that the SkyView cache is *not* the same as the ``fcm_data`` folder in which fcmaker stores the downloaded FITS files.

For example, if you want to boost the aesthetic appeal of your finding charts with the 
exquisite Computer Modern Bright sans-serif font (why wouldn't you?), try::

   python -m fcmaker --systemtex

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
   * **XSHOOTER**: ``DSS2 Red``
   
.. note::
   
   From an operational perspective, when using Skyview images, I would strongly recommend 
   to only use ``DSS2 Red`` or ``SDSSr`` for optical instruments, and ``2MASS-J``, 
   ``2MASS-H`` or ``2MASS-K`` for IR instruments.

