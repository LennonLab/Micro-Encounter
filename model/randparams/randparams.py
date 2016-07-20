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

    width  = 10 #randint(5, 10)
    height = 10 #randint(5, 10)

    num_envgrads = 1
    for i in range(num_envgrads):
        x = np.random.uniform(0, width)
        y = np.random.uniform(0, height)
        envgrads.append([x, y])

    m = 0.0
    speciation = 0.0
    alpha = 0.99
    r = choice([1, 2, 3, 4])#, 5, 6, 7, 8, 9, 10])#, 12, 14, 16, 20, 24, 28, 32, 40, 50, 70, 90]) # resource particles flowing in per time step

    gmax = np.random.uniform(0.1, 0.5) # max specific growth rate
    maintmax = np.random.uniform(0.025, 0.05)

    dmax = np.random.uniform(0.01, 0.1)  # max dispersal probability
    pmax = np.random.uniform(0.001, 0.01)  # max probability of going active
    mmax = randint(30, 50)  # max maintenance factor
    std = np.random.uniform(0.1, 0.4)
    maxgen = randint(1, 10)

    # TO EXPLORE A SINGLE SET OF VALUES FOR MODEL PARAMETERS

    plist = [width, height, alpha, speciation, seedCom, m, r, gmax, maintmax, dmax, envgrads, rates, pmax, mmax, std, maxgen]
    return plist
