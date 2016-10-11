from __future__ import division
from random import randint
import numpy as np
import sys

def get_rand_params(extremes):
    """ Get random model parameter values. Others are chosen in bide.py """

    seed = 100 # size of starting community
    dim = randint(1000, 1000)
    width  = int(dim) # in microns
    height = int(dim) # in microns
    length = int(dim) # in microns

    m = 0.0001 #np.random.uniform(0.01, 0.001)
    r = 0.01 #np.random.uniform(0.01, 0.1) # resource particles flowing in per time step
    std = np.random.uniform(1, 1)

    gmax = 0.5 #np.random.uniform(0.1, 1.0) # max specific growth rate
    dmax = 0.5 #np.random.uniform(0.1, 1.0)  # max dispersal rate
    pmax = 0.5 #np.random.uniform(0.1, 1.0)  # max probability of going active

    mmax = 4 #np.random.uniform(1, 10)
    mfact = 20 #randint(1, 30) # max maintenance factor

    if extremes == True:

        x = np.random.binomial(1, 0.5)
        if x == 0: # low dormant capacity
            pmax = 0.1
            mfact = 10
            plist = [width, height, length, seed, m, r, gmax, mmax, dmax, pmax, mfact, std]

        elif x == 1: # high dormant capacity
            pmax = 0.001
            mfact = 1000
            plist = [width, height, length, seed, m, r, gmax, mmax, dmax, pmax, mfact, std]

        return plist


    else:
        plist = [width, height, length, seed, m, r, gmax, mmax, dmax, pmax, mfact, std]
        return plist
