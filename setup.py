#!/usr/bin/env python
from setuptools import setup

setup(name='backstop_history',
      packages=["backstop_history"],
      use_scm_version=True,
      setup_requires=['setuptools_scm', 'setuptools_scm_git_archive'],
      description='Tools to assemble backstop command histories for Thermal Models',
      url='http://github.com/acisops/backstop_history',
      author='Gregg Germain',
      author_email='ggermain@cfa.harvard.edu',
      include_package_data=True,
      zip_safe=False)
