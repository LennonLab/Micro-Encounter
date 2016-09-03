from __future__ import division
from random import randint, seed
import numpy as np

def get_rand_params():
    """ Get random model parameter values. Others are chosen in bide.py """

    seed()
    seedCom = 100 # size of starting community
    width  = 10 #randint(5, 10)
    height = 10 #randint(5, 10)

    m = 0 #np.random.uniform(0.01, 0.1)
    r = randint(1, 1) # resource particles flowing in per time step
    std = np.random.uniform(0.4, 0.4)

    gmax = np.random.uniform(50, 50) # max specific growth rate
    dmax = np.random.uniform(0.01, 0.01)  # max dispersal probability
    pmax = np.random.uniform(0.01, 0.01)  # max probability of going active

    maintmax = np.random.uniform(10, 100)
    mmax = randint(10, 10)  # max maintenance factor

    plist = [width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std]
    return plist
