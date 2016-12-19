#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 00:31:09 2016

@author: gionuno
"""

import cython
from libcpp.vector cimport vector

cdef extern from "liblinprog.hpp":
    cdef double pdim(vector[vector[double]] &,vector[double] &,vector[double] &,vector[double] &);
    
def interior_method(a,b,c,eq):
    return pdim(a,b,c,eq);