
Installation
============

fcmaker is available on testpypi, which makes its installation easier than ever. 
In a terminal, type:
::

   pip install fcmaker

And that should take care of things.

The most recent release of fcmaker is also available for download from its `Github repository <https://github.com/fpavogt/fcmaker/releases/latest/>`_. 
Interested users can fork the fcmaker repository if they want to get access to the 
latest updates not yet released. *Push requests* for bug fixes and new features are 
welcome and will be examined in detail. 
      
Requirements
------------
fcmaker is written in Python 3.6. The following packages are required for it to work 
properly:

* numpy (1.13.1 or or above)
* scipy (0.19.0 or above)
* matplotlib (2.0.2 or above)
* astropy (2.0.1 or above)
* astroquery (0.3.4 or above)
* aplpy (1.1.1 or above)
* p2api (0.9 or above)

Optional (but strongly recommended): 

* Montage and montage-wrapper (0.9.9 or above)
* Proper LateX installation

The Montage package is required to rotate the finding charts North, even if the underlying
FITS file isn't. A system-wide LaTeX installation allows for a slick(-er) looking font for
the various labels.

Testing the installation
------------------------

In a terminal shell, try to access the basic help of fcmaker::
 
   python -m fcmaker --help
 
If that works, chances are, you will probably be fine.

.. _troubleshooting:

Troubleshooting
---------------

1. 
If you encounter errors when running fcmaker, first ensure that your numpy, scipy, 
matplotlib, aplpy, and astropy packages are up-to-date, then try again. 
  
2. 
If you still encounter errors after doing that, check the :ref:`faq`.
  
3. 
Check if this is a known issue: https://github.com/fpavogt/fcmaker/issues
  
4. 
If you still can't figure out what's wrong, please `submit a new issue on the Github 
repository of the project <https://github.com/fpavogt/fcmaker/issues>`_. Provide as much 
detail as possible (error message, minimal example able to reproduce the error, 
operating system, Python version, etc ...).

.. note::
   Submitting a Github issue is the best way for you to get help rapidly, for us to keep 
   track of the problems that need solving, and for future users to see what changes have 
   been made over time (and look for existing solutions to their problem which may be the 
   same as yours). Submitting a new issue on Github is rapid and easy, but if you are 
   really against doing it (why would you ?), you can always email 
   frederic.vogt@alumni.anu.edu.au for help. 

 