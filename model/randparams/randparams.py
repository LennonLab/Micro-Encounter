from __future__ import division
from random import choice, randint
import numpy as np
import sys
import os

def get_rand_params(fixed):
    """ Get random model parameter values. Others are chosen in bide.py """

    envgrads = []
    seedCom = 10 # size of starting community
    rates = []

    rates = np.array([0.1])

    width  = 20
    height = 20

    num_envgrads = 1
    for i in range(num_envgrads):
        x = np.random.uniform(0, width)
        y = np.random.uniform(0, height)
        envgrads.append([x, y])

    m = 0.0
    speciation = 0.0
    alpha = 0.99
    r = 6 #choice([1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]) # resource particles flowing in per time step

    gmax = 0.2 #np.random.uniform(0.01, 0.1) # max specific growth rate
    maintmax = 0.001 #np.random.uniform(0.001, 0.005)

    dmax = 0.03 #np.random.uniform(0.001, 0.01)  # max dispersal probability
    pmax = 0.2 #np.random.uniform(0.01, 0.1)  # max probability of going active
    mmax = 20 #randint(20, 40)  # max maintenance factor
    std = 0.4 #np.random.uniform(0.01, 0.4)

    # TO EXPLORE A SINGLE SET OF VALUES FOR MODEL PARAMETERS

    plist = [width, height, alpha, speciation, seedCom, m, r, gmax, maintmax, dmax, envgrads, rates, pmax, mmax, std]
    return plist
