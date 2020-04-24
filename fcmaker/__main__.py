# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0311
'''
 fcmaker: a Python module to automatically create finding charts for ESO OBs in p2.\n
 Copyright (C) 2017-2019,  F.P.A. Vogt
 --- oOo ---

 This file contains the master fcmaker routines. See the dedicated website for more info:
 http://fpavogt.github.io/fcmaker
 Created October 2017, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au

 --- oOo ---

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
# --------------------------------------------------------------------------------------------------

# Import generic python packages
import os
import argparse
from datetime import datetime

import yaml

# Import the fcmaker submodules
from .fcmaker_version import __version__
from . import fcmaker as fcm

# Use argparse to make fcmaker user friendly -------------------------------------------------------
parser = argparse.ArgumentParser(description='''Creates ESO-compliant finding charts from
                                                OBs on p2 (or locally from a parameter file).
                                                If no parameter file is passed, connect
                                                to p2 and select an obId manually. ''',
                                 epilog='Full documentation: %s \n \n \
                                         Feedback, questions, comments: \
                                         frederic.vogt@alumni.anu.edu.au \n' %
                                         ('http://fpavogt.github.io/fcmaker'),
                                 formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--version', action='version', version=('fcmaker %s'%__version__))

parser.add_argument('--p2uid', action='store', metavar='your_p2_username', default=None, help='p2 user ID')

parser.add_argument('--demo', action='store_true', help='Connect to the demo p2 server')

parser.add_argument('--obid', action='store', metavar='12345678', default=None,
                    help='Observing Block ID on p2')

parser.add_argument('-f', action='store', metavar='filename', nargs=1,
                    type=open, help='parameter filename')

parser.add_argument('--obsdate', action='store', metavar='%Y-%m-%d %H:%M:%S', nargs='+',
                    default=[datetime.strftime(datetime.utcnow(), '%Y-%m-%d %H:%M:%S') + ' UTC'],
                    help='Date of the observations (for targets with proper motions)')

parser.add_argument('--bk-image', action='store', metavar='some_image.fits', nargs='+',
                    default=None, help='filename for the background image')

parser.add_argument('--bk-lam', action='store', metavar='inst + lambda', nargs='+',
                    default=None, help='Wavelength of the background image')

parser.add_argument('--do-pdf', action='store_true',
                    help='save a pdf version of the chart (in addition to the jpg)')

parser.add_argument('--do-png', action='store_true',
                    help='save a png version of the chart (in addition to the jpg)')

parser.add_argument('--do-parang', action='store_true',
                    help='If a parallactic angle is requested, print the instrument field-of-view')

parser.add_argument('--data-loc', action='store', metavar='folder/name/',
                    default=os.path.join('.', 'fcm_data'), help='Location to store the data')

parser.add_argument('--plot-loc', action='store', metavar='folder/name',
                    default=os.path.join('.', 'fcm_plots'), help='Location to store the charts')

parser.add_argument('--no-upload', action='store_true',
                    help='disable the upload of finding charts to p2')

parser.add_argument('-l', '--local', action='store_true', help='feed a local OB description')

parser.add_argument('--montage', action='store_true',
                    help='enable the use of Montage to rotate the charts')

parser.add_argument('--systemtex', action='store_true',
                    help='disable the use of the system-wide LaTeX')

parser.add_argument('--clear-SkyView-cache', action='store_true', help='clear the SkyView cache')


def main():
   ''' The main function, runs when fcmaker is started with the high-level entry point. '''

   # What did the user type in ?
   args = parser.parse_args()

   # Ok, do I need to deal with a local OB file ?
   if args.local:
   # Yes!

      # Did I get a parameter file ?
      if args.f is None: # No !
         raise Exception('Ouch! To use the local mode, you must specify a suitable '+
                         'parameter file with "-f filename".')

      fcm.make_fc_local(args.f[0],
                        do_pdf=args.do_pdf,
                        do_png=args.do_png,
                        systemtex=args.systemtex,
                        montage=args.montage,
                        clear_SkyView_cache=args.clear_SkyView_cache,
                        obsdate=' '.join(args.obsdate),
                        do_parang=args.do_parang,
                        )

   else:
   # Ok, we'll look on P2 for the info ...

      # Did I get a parameter file ?
      if args.f is not None: # Yes !

         # Load the parameter file
         inpars = yaml.load(args.f[0], Loader=yaml.SafeLoader)

         # If I need to connect to p2, extract the user ID and obIDs.
         p2uid = inpars['p2uid']
         pswd = inpars['pswd']
         obids = inpars['obids']
         data_loc = inpars['data_loc']
         plot_loc = inpars['plot_loc']
         bk_images = inpars['bk_images']
         bk_lams = inpars['bk_lams']

      else: # ok, a manual run it is !

         # Get the p2 login info
         p2uid = args.p2uid
         pswd = None

         # What OB should we create a finding hart for ?
         if args.obid is not None:
            obids = [int(args.obid)]
         else:
            obids = []

         # Just use the default bk_image in manual mode
         if args.bk_image is not None:
            bk_images = [' '.join(args.bk_image)]
         else:
            bk_images = [args.bk_image]

         if args.bk_lam is not None:
            bk_lams = [' '.join(args.bk_lam)]
         else:
            bk_lams = [args.bk_lam]

         data_loc = args.data_loc
         plot_loc = args.plot_loc

      # Launch the main fcmaker routine
      fcm.make_fc(p2uid=p2uid, pswd=pswd, demo=args.demo, obids=obids,
                  bk_images=bk_images, bk_lams=bk_lams,
                  data_loc=data_loc, plot_loc=plot_loc,
                  do_pdf=args.do_pdf, do_png=args.do_png,
                  no_upload=args.no_upload, systemtex=args.systemtex,
                  montage=args.montage, clear_SkyView_cache=args.clear_SkyView_cache,
                  obsdate=' '.join(args.obsdate), do_parang=args.do_parang)

# Start of the interactive part --------------------------------------------------------------------
if __name__ == "__main__":

   main()

# ----------------- End of the World as we know it -------------------------------------------------
# --------------------------------------------------------------------------------------------------
