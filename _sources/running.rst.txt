.. _examples:
.. highlight:: python

Usage instructions
==================

There are two main ways to use fcmaker. You can either execute the module as a script, i.e.::

   python -m fcmaker
   
or you can import the module and execute it within a script, i.e.::

   >>> import fcmaker
   >>> fcmaker.make_fc( ... ) # or fcmaker.make_fc_local ( ... )
   
When running fcmaker as a script, any argument you feed it is pretty much sent to the 
underlying functions ``make_fc()`` and ``make_fc_local()``. For simplicity, this page
only discusses how to run the entire module as a script, which ought to be slightly more 
friendly to users not (yet!) familiar with Python. Still, the hope is that after reading 
this page, the use of the functions ``make_fc()`` and ``make_fc_local()`` shouldn't be
too mysterious. Also, don't forget the built-in help::

   >>> import fcmaker
   >>> help(fcmaker.make_fc)
   

1: creation of a finding chart for a single OB on p2
-----------------------------------------------------------

In its most basic mode, fcmaker connects to `p2 <http://www.eso.org/p2>`_ and creates the 
finding chart for a given OB Id. To do so, simply type in a  terminal::

   python -m fcmaker
   
You will be prompted for your `p2 <http://www.eso.org/p2>`_ username, password, and the 
ID of the observing block to process. fcmaker will create two folders ``fcm_data`` and 
``fcm_plots`` at your current location, where it will store the background image and the 
finding chart. Once the finding chart has been generated, it will either be attached to 
the OB directly (if the OB contains no pre-existing chart), or you will be asked to select
an upload slot (occupied, or not). **Note that an OB cannot contain two finding charts with 
the same name!**

If you get fed up with typing your p2 username and/or the OB id repeatedly, you can 
specify them directly when calling fcmaker::

   python -m fcmaker --p2uid jdoe --obid 12345678 


2: creation of a finding chart for a single OB on p2demo
--------------------------------------------------------------

If you simply want to experiment with p2 and/or fcmaker, you can use the `p2demo <http://www.eso.org/p2>`_
server to create mock OBs. In that case, simply use the ``--demo`` flag to tell fcmaker that 
you want to connect to the tutorial account on the demo server::

   python -m fcmaker --demo

   
You will only be prompted for the ID of the observing block to process.

3: creation of finding charts for many OBs on p2 (with upload)
------------------------------------------------------------------------

If you have a lot of finding charts to create, fcmaker allows for batch processing.
First, create a text file :download:`p2_2_fcm.params.txt <./examples/p2_2_fcm.params.txt>` 
(the actual filename is flexible) with the following structure:

.. literalinclude:: examples/p2_2_fcm.params.txt

Then, feed it to fcmaker with the ``-f`` flag::

   python -m fcmaker -f p2_2_fcm.params.txt

In doing so, fcmaker will connect to `p2 <http://www.eso.org/p2>`_ and process all the 
OBs listed. Note that for each finding chart, you will still need to manually specify 
whether you want to upload it to `p2 <http://www.eso.org/p2>`_ (or not) if the OB already 
contains at least 1 finding chart. 

4: creation of finding charts for many OBs on p2 (no upload)
-----------------------------------------------------------------------

You can fully automate the creation of many finding charts if you include the 
``--no-upload`` flag. In that case, fcmaker will never attempt to upload anything to 
`p2 <http://www.eso.org/p2>`_, and only save the image files locally::

   python -m fcmaker --no-upload -f p2_2_fcm.params.txt


6: specifying observing dates and times
---------------------------------------

For targets with large proper motions, one can provide the expected year, (month, 
day, ...!) of the observation to create an accurate finding chart, using the ``--obsdate`` 
flag. This also works for OBs that have Ephemeris files. For example::

   python -m fcmaker -f --obsdate 2018
   python -m fcmaker -f --obsdate 2018-05
   python -m fcmaker -f --obsdate 2018-05-17 14:34:57 UTC
   
Finding charts for moving targets get automatically tagged with the symbol :math:`\leadsto`.
   
   
7: creation of finding charts locally (without p2)
-------------------------------------------------------

fcmaker can also create finding charts *locally*, without the need to have an OB present 
on `p2 <http://www.eso.org/p2>`_ first (fcmaker will still require an internet connection, 
though!). To do so, one first needs to specify the basic OB parameters in a text file  
(in essence, a stripped-down version of a full OBX file). Here are
templates for all supported instruments:
   
   - :download:`local_2_fcm.muse.txt  <./examples/local_2_fcm.muse.txt>`
   - :download:`local_2_fcm.hawki.txt <./examples/local_2_fcm.hawki.txt>`
   - :download:`local_2_fcm.xshooter.txt <./examples/local_2_fcm.xshooter.txt>`
   - :download:`local_2_fcm.espresso.txt <./examples/local_2_fcm.espresso.txt>`
   
The file is then fed to fcmaker with the ``-f`` flag, together with the ``--local`` flag to 
indicate that it is a *local run*::

   python -m fcmaker --local -f local_2_fcm.muse.txt
   python -m fcmaker --local -f local_2_fcm.hawki.txt
   python -m fcmaker --local -f local_2_fcm.xshooter.txt
   python -m fcmaker --local -f local_2_fcm.espresso.txt
   
fcmaker will create the associated finding chart, store it where specified, and exit.


8: the background images
------------------------

With the exception of MUSE NFM and ESPRESSO finding charts (see :ref:`gaia-images`), fcmaker relies on 
``astroquery.skyview`` to download background images for the finding charts (unless a local 
FITS file is provided by the user, see :ref:`local-FITS`). To display the full list of 
surveys available, type in a Python shell::
 
   from astroquery.skyview import SkyView
   SkyView.survey_dict['overlay_blue']

The default background survey images for the instruments supported by fcmaker are as follows:
   
   * **MUSE WFM**: ``DSS2 Red``
   * **HAWKI**: ``2MASS-J``, ``2MASS-H`` and ``2MASS-K`` for filters ``J``, ``H``, ``K``, and ``2MASS-H`` for all other filters
   * **XSHOOTER**: ``DSS2 Red``
   
.. note::
   
   From an operational perspective, when using Skyview images, I would strongly recommend 
   to only use ``DSS2 Red`` or ``SDSSr`` for optical instruments, and ``2MASS-J``, 
   ``2MASS-H`` or ``2MASS-K`` for IR instruments.


.. _local-FITS: 

Using local FITS files
++++++++++++++++++++++

Users interested in providing their own background image can do so (for the left-hand-side plot)
via the ``--bk-image`` and ``--bk-lam`` parameters. The former specifies the filename, 
assuming that the FITS file is placed in the ``--data-loc`` directory (default: ``./fcm_data/``). 
The latter allows the user to specify the text to be displayed on the chart. To meet ESO's 
requirements, this text should specify the wavelength of the image::

   python -m fcmaker --p2uid joe --obid 12345678 --bk-image my-file.fits --bk-lam 450-550nm (HST)

See the Gallery (:ref:`examples-NFM`) for an illustration of this approach.


.. _gaia-images:

Mock Gaia images for MUSE NFM
+++++++++++++++++++++++++++++

``DSS2 Red`` background images are not well suited for MUSE NFM finding charts, given the small
field-of-view of 7.5x7.5 square arcsec of this mode, or for ESPRESSO finding charts, given the 
brightness of the targets. To circumvent this issue, fcmaker can create pseudo sky images 
from the Gaia catalogue, via the ``fcmaker_plots.make_gaia_image()`` function. By default, 
the image is created with a pixel size of 25 mas for MUSE NFM (the pixel size) and 50 mas for ESPRESSO. 
Each star is plotted as a 2D gaussian with a FWHM of 80 mas for MUSE NFM (typical to the image 
quality achieved in normal operations) or 0.6 arcsec for ESPRESSO (typical V-band seeing), and 
scaled in intensity as a function of its Gaia flux. The proper motion of each star is 
propagated up to the OB observing date (requested by the user with the ``--obsdate`` parameter; default = chart 
creation time). The resulting image is saved as a fully-fledged FITS file at the 
``--data-loc`` location (default: ``./fcm_data/``).  

See the Gallery (:ref:`examples-NFM`) for an illustration of the benefit of these mock 
Gaia images.


.. _flags:

The fcmaker flags
-----------------

Here are all the flags that allow to fine-tune the way fcmaker works. They are:
 
 * ``--help,-h`` : prints the basic help
 * ``--version`` : prints the fcmaker version
 * ``--p2uid username``: p2 user ID
 * ``--obid 12345678``: Observing Block ID on p2
 * ``--demo``: Connects to the p2demo site rather than the production site. Overrides p2uid and pswd.
 * ``-f filename``: specify a parameter filename
 * ``--obsdate YYYY-MM-DD hh:mm:ss TZ`` : allows to specify the observing date (and time, and time zone), in case of large proper motions or moving target. Assuming UTC by default.
 * ``--bk-image filename.fits``: a FITS filename for the background image. Must be stored in ``--data-loc``.
 * ``--bk-lam some text``: specifies the wavelength of the background image.
 * ``--do-pdf`` : tells fcmaker to save a .pdf version of the finding chart (in addition to the jpg one).
 * ``--do-png`` : tells fcmaker to save a .png version of the finding chart (in addition to the jpg one).
 * ``--do-parang``: forces the drawing of the instrument field-of-view when a parallactic angle is set. This will be accurate, but *very* time-dependant!
 * ``--data-loc path/to/folder``: relative path to store the downloaded data (default: ./fcm_data)
 * ``--plot-loc path/to/folder``: relative path to store the finding charts (default: ./fcm_plots)
 * ``--no-upload`` : tells fcmaker never to try to upload the finding chart(s) to p2. 
 * ``-l`` or ``--local``: flags the parameter file as a local OB definition
 * ``--montage`` : will use Montage, to force the finding charts to be rotated North (soon to be retired).
 * ``--systemtex`` : will use the system-wide LaTeX, rather than the matplotlib one. Finding charts will be prettier. [encouraged!]
 * ``--clear-SkyView-cache`` : does as it says. Note that the SkyView cache is *not* the same as the ``fcm_data`` folder in which fcmaker stores the downloaded FITS files.


For example, if you want to boost the aesthetic appeal of your finding charts with the 
exquisite Computer Modern Bright sans-serif font (why wouldn't you?), try::

   python -m fcmaker --systemtex

