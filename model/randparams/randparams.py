from __future__ import division
from random import randint, seed
import numpy as np

def get_rand_params():
    """ Get random model parameter values. Others are chosen in bide.py """

    seed()
    seedCom = 100 # size of starting community
    width  = 10 #randint(5, 10)
    height = 10 #randint(5, 10)

    m = np.random.uniform(0.001, 0.001)
    r = randint(2, 2) # resource particles flowing in per time step
    std = np.random.uniform(0.4, 0.4)

    gmax = np.random.uniform(0.1, 0.1) # max specific growth rate
    dmax = np.random.uniform(0.5, 0.5)  # max dispersal rate
    pmax = np.random.uniform(0.01, 0.01)  # max probability of going active

    maintmax = np.random.uniform(10, 10)
    mmax = randint(20, 20)  # max maintenance factor

    plist = [width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std]
    return plist
