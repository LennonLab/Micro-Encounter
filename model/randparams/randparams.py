from __future__ import division
from random import randint
import numpy as np
import sys

def get_rand_params(extremes):
    """ Get random model parameter values. Others are chosen in bide.py """

    seed = 100 # size of starting community
    dim = randint(43200, 43200)
    width  = int(dim) # in microns
    height = int(dim) # in microns
    length = int(dim) # in microns

    m = np.random.uniform(0.001, 0.01)
    r = np.random.uniform(0.01, 0.1) # resource particles flowing in per time step
    std = np.random.uniform(1, 1)

    gmax = np.random.uniform(0.1, 1.0) # max specific growth rate
    dmax = np.random.uniform(0.1, 1.0)  # max dispersal rate
    pmax = np.random.uniform(0.001, 0.01)  # max probability of going active

    mmax = np.random.uniform(1, 10)
    mfact = randint(10, 40) # max maintenance factor

    plist = [width, height, length, seed, m, r, gmax, mmax, dmax, pmax, mfact, std]
    return plist
