from setuptools import setup, find_packages  # Always prefer setuptools over distutils
import os

# Run the version file
v = open(os.path.join('.','fcmaker','fcmaker_version.py'))
version = [l.split("'")[1] for l in v.readlines() if '__version__' in l][0]

setup(
   name='fcmaker',
   version=version,
   author='F.P.A. Vogt',
   author_email='frederic.vogt@alumni.anu.edu.au',
   packages=['fcmaker',],
   url='http://fpavogt.github.io/fcmaker/',
   download_url='https://github.com/fpavogt/fcmaker/archive/master.zip',
   license='GNU General Public License',
   description='Python module to automatically create ESO-compliant finding charts in p2.',
   long_description=open('README').read(),
   python_requires='>=3',
   install_requires=[
      "numpy >= 1.16.2",
      "scipy >= 1.2.1",
      "matplotlib >= 3.0.3",
      "astropy >=3.1.2",
      "aplpy >=2.0.3",
      "astroquery >= 0.3.9",
      "p2api >= 0.94",
      "PyYAML >=5.1",
      "pytz >= 2018.9",
      "astroplan >=0.4",
      "pillow >=5.4.1",
      #"montage-wrapper >= 0.9.9",
   ],
   
   entry_points={'console_scripts': ['fcmaker=fcmaker.__main__:main']},

   classifiers=[
   # How mature is this project? Common values are
   #   3 - Alpha
   #   4 - Beta
   #   5 - Production/Stable
   'Development Status :: 4 - Beta',

   # Indicate who your project is intended for
   'Intended Audience :: Science/Research',
   'Topic :: Scientific/Engineering :: Astronomy',

   # Pick your license as you wish (should match "license" above)
   'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
    
   # Specify the Python versions you support here. In particular, ensure
   # that you indicate whether you support Python 2, Python 3 or both.
   'Programming Language :: Python :: 3.7',
   ],
    
   include_package_data=True, # So that non .py files make it onto pypi, and then back !
   #package_data={
   #     'example_files':['example_files/*'],
   #     'docs':['../docs/build']
   #    }
)