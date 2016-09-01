from __future__ import division
from random import randint, seed
import numpy as np

def get_rand_params():
    """ Get random model parameter values. Others are chosen in bide.py """

    seed()
    seedCom = 10 # size of starting community
    width  = 10 #randint(5, 10)
    height = 10 #randint(5, 10)

    m = 0 #np.random.uniform(0.01, 0.1)
    r = randint(1, 10) # resource particles flowing in per time step

    gmax = np.random.uniform(0.1, 0.1) # max specific growth rate
    maintmax = np.random.uniform(0.01, 0.01)

    dmax = np.random.uniform(0.01, 0.1)  # max dispersal probability
    pmax = np.random.uniform(0.1, 0.1)  # max probability of going active
    mmax = randint(1, 1)  # max maintenance factor
    std = np.random.uniform(0.04, 0.4)

    # TO EXPLORE A SINGLE SET OF VALUES FOR MODEL PARAMETERS

    plist = [width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std]
    return plist
