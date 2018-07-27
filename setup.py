from setuptools import setup, find_packages  # Always prefer setuptools over distutils

setup(
   name='fcmaker',
   version='0.3.8',
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
      "numpy >= 1.12",
      "scipy >= 0.19.0",
      "matplotlib >= 2.0.2",
      "astropy >=3.0",
      "aplpy >=1.1.1",
      "astroquery >= 0.3.4",
      "p2api >= 0.92",
      "PyYAML >=3.12",
      "pytz >= 2018",
      "astroplan",
      "pillow >=4.2.1",
      #"montage-wrapper >= 0.9.9",
   ],
    
   classifiers=[
   # How mature is this project? Common values are
   #   3 - Alpha
   #   4 - Beta
   #   5 - Production/Stable
   'Development Status :: 3 - Alpha',

   # Indicate who your project is intended for
   'Intended Audience :: Science/Research',
   'Topic :: Scientific/Engineering :: Astronomy',

   # Pick your license as you wish (should match "license" above)
   'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
    
   # Specify the Python versions you support here. In particular, ensure
   # that you indicate whether you support Python 2, Python 3 or both.
   'Programming Language :: Python :: 3.6',
   ],
    
   include_package_data=True, # So that non .py files make it onto pypi, and then back !
   #package_data={
   #     'example_files':['example_files/*'],
   #     'docs':['../docs/build']
   #    }
)