#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 17 19:48:00 2016

@author: gionuno
"""

import matplotlib.pyplot as plt;
import matplotlib.image  as img;
import numpy as np;

from liblinprog import interior_method;

BINS = 4;
hist = np.load("hist_"+str(BINS)+".npy");

n = BINS**3;
c = np.zeros((n*n));
for i1 in range(BINS):
    for i2 in range(BINS):
        for i3 in range(BINS):
            for j1 in range(BINS):
                for j2 in range(BINS):
                    for j3 in range(BINS):
                        c[n*(BINS*BINS*i1+BINS*i2+i3)+BINS*BINS*j1+BINS*j2+j3] = 1.0*(abs(i1-j1)+abs(i2-j2)+abs(i3-j3));

O1 = np.ones((1,n));                         
O2 = np.eye(n);

A = np.concatenate((np.kron(O2,O1),np.kron(O1,O2),np.ones((1,n*n))),axis=0);

b = np.zeros(2*n+1);
b[-1] = 1.0;

e = -np.ones(2*n+1);
e[-1] = 0.0;
       
dist_e = np.zeros((len(hist),len(hist)));
dist_h = np.zeros((len(hist),len(hist)));
for i in range(len(hist)):
    for j in range(i+1,len(hist)):
        print i,j;
        b[:n] = np.copy(hist[i,:]);
        b[n:-1] = np.copy(hist[j,:]);
        dist_e[i,j] = dist_e[j,i] = interior_method(A,b,c,e);
        dist_h[i,j] = dist_h[j,i] = np.sqrt(1.0-np.sum(np.sqrt(hist[i,:]*hist[j,:])));

plt.imshow(dist_h);
plt.show();
plt.imshow(dist_e);
plt.show();

np.save("dist_e_"+str(BINS),dist_e);
np.save("dist_h_"+str(BINS),dist_h);
