#!/usr/bin/env python
from setuptools import setup
from backstop_history import __version__


setup(name='backstop_history',
      packages=["backstop_history"],
      version=__version__,
      description='Tools to assemble backstop command histories for Thermal Models',
      url='http://github.com/acisops/backstop_history',
      download_url=url,
      author='Gregg Germain',
      author_email='ggermain@cfa.harvard.edu',
      include_package_data=True,
      classifiers=[
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
      ],
      zip_safe=False)
