#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 17 19:48:00 2016

@author: gionuno
"""

import matplotlib.pyplot as plt;
import matplotlib.image  as img;
import numpy as np;

import io;
import re;
import os;

import glob;
import string;

names = glob.glob("images/*.jpg");

def get_hist(I,B):
    a = np.zeros([B*B*B]);
    for i in range(I.shape[0]):
        for j in range(I.shape[1]):
            r = int((I[i,j,0]*B)/256);
            g = int((I[i,j,1]*B)/256);
            b = int((I[i,j,2]*B)/256);
            a[B*B*r+B*g+b] += 1.0;
    return a / (I.shape[0]*I.shape[1]);

BINS = 4;

H = [];
S = [];
for name in names:
    print name;
    S.append(name);
    H.append(get_hist(img.imread(name),BINS));

S = np.array(S);
H = np.array(H);
        
print 'Saving';
np.save("hist_"+str(BINS),H);
np.save("hist_names_"+str(BINS),S);
print 'Done';