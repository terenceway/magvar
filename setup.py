#!/usr/bin/env python3

from distutils.core import setup, Extension

magvar_module = Extension('magvar',
			  sources=['av_mag.c', 'magvarmodule.c', 'wmm.c'])

setup(name='magvar',
      version='1.0.0',
      description='Magnetic variation (declination) at a given date, latitude, longitude, and elevation',
      ext_modules = [magvar_module],
      author='Terence Way',
      author_email='terry@wayforward.net')
