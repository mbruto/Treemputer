#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-

import setuptools
from distutils.core import Extension

C_extension = Extension(

 name = 'Library_iaf',
 sources = ["libc/Library_iaf.c"],
 library_dirs = ["libc"],
 include_dirs = ["libc"]
 
)



setuptools.setup(name='Iterative_additive_filling',
 version='1.0',
 description='Impute missing distances in (nearly) additive sparse matrices ',
 author='Maxime Bruto',
 author_email='bruto.maxime@gmail.com',
 scripts=['Iterative_additive_filling.py'],
 install_requires=['argparse', 'numpy'],
 packages=setuptools.find_packages(),
 include_package_data=True,
 ext_modules=[C_extension]
)
