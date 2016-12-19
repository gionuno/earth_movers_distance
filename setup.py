#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 00:38:51 2016

@author: gionuno
"""

from distutils.core import setup;
from distutils.extension import Extension;
from Cython.Build import cythonize;

e_modules = cythonize([Extension("liblinprog",["liblinprog.pyx"],libraries=["armadillo"],language="c++")]);
setup(name="liblinprog",ext_modules = e_modules)