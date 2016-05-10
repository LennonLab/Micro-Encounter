from __future__ import division
from random import choice, randint
import numpy as np
import sys
import os

def get_rand_params(fixed):
    """ Get random model parameter values. Others are chosen in bide.py """

    envgrads = []
    seedCom = 100 # size of starting community
    rates = []

    if fixed is True:

        #rates = np.array([1.0, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001])
        #rates = [round(10**np.random.uniform(-2, 0.0), 4)]
        rates = np.array([0.1])

        width  = 20
        height = 20

        num_envgrads = 2
        for i in range(num_envgrads):
            x = np.random.uniform(0, width)
            y = np.random.uniform(0, height)
            envgrads.append([x, y])

        m = 0.5
        speciation = 0.01
        alpha = 0.99
        r = 10

        gmax = 0.2 # max specific growth rate
        maintmax = 0.01

        rmax = 200 # max resource particle size
        dmax = 0.01 # max dispersal probability
        pmax = 0.01 # max probability of going active
        mmax = 10 # max maintenance factor


    elif fixed is False:

        rates = np.array([1.0, 0.7, 0.4, 0.1, 0.07, 0.04, 0.01, 0.007])

        #width = randint(6, 10)
        width = 100

        #height = randint(10, 100)
        height = 100

        alpha = np.random.uniform(0.95, 0.99)
        speciation = np.random.uniform(0.01, 0.1)
        m = np.random.uniform(0.0001, 0.001) # m = probability of immigration
        m = 0.0

        r = randint(1, 10) #resource particles flowing in per time step
        rmax = randint(10, 100) # maximum resource particle size
        rmax = 100

        num_envgrads = randint(1, 2)
        for i in range(num_envgrads):
            x = np.random.uniform(0, width)
            y = np.random.uniform(0, height)
            envgrads.append([x, y])

        gmax = np.random.uniform(0.05, 0.5)
        dmax = np.random.uniform(0.01, 0.1) # probability of dispersing in a given time step
        maintmax = np.random.uniform(0.001*gmax, 0.01*gmax) # maximum metabolic maintanence cost
        pmax = np.random.uniform(0.5, 0.9) # max probability of going active
        mmax = np.random.uniform(1, 2) # max maintenance 'factor'


        # TO EXPLORE A SINGLE SET OF VALUES FOR MODEL PARAMETERS

    return [width, height, alpha, speciation, \
            seedCom, m, r, rmax, gmax, maintmax, dmax, \
            envgrads, rates, pmax, mmax]
