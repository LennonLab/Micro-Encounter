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

        rates = np.array([1.0, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001])

        motion = 'white_noise'
        #motion = 'brown_noise'

        width  = 20
        height = 20

        num_envgrads = 2
        for i in range(num_envgrads):
            x = np.random.uniform(0, width)
            y = np.random.uniform(0, height)
            envgrads.append([x, y])

        nNi = 1 # max number of Nitrogen types
        nP = 1 # max number of Phosphorus types
        nC = 1 # max number of Carbon types

        m = 0.0
        speciation = 0.01
        alpha = 0.98
        r = 10

        gmax = 0.1 # max specific growth rate
        maintmax = 0.1*gmax

        rmax = 100 # max resource particle size
        dmax = 0.01 # max dispersal probability
        pmax = 0.01 # max probability of going active
        mmax = 100 # max maintenance factor


    elif fixed is False:

        #motion = choice(['fluid', 'white_noise']) # 'fluid', 'unidirectional'
        motion = 'fluid'
        #motion = 'brown_noise'

        if motion == 'white_noise':
            rates = np.array([0.00001])
        else:
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

        nNi = randint(1, 10) # max number of Nitrogen types
        nP = randint(1, 10) # max number of Phosphorus types
        nC = randint(1, 10) # max number of Carbon types

        num_envgrads = randint(1, 10)
        for i in range(num_envgrads):
            x = np.random.uniform(0, width)
            y = np.random.uniform(0, height)
            envgrads.append([x, y])

        #tp_max = np.random.uniform(0.001, 0.1) # maximum probability of transitioning at random from dormant to active (i.e., Scout hypothesis)

        gmax = np.random.uniform(0.05, 0.5)
        dmax = np.random.uniform(0.01, 0.1) # probability of dispersing in a given time step
        maintmax = np.random.uniform(0.001*gmax, 0.01*gmax) # maximum metabolic maintanence cost
        pmax = np.random.uniform(0.5, 0.9) # max probability of going active
        mmax = np.random.uniform(1, 2) # max maintenance 'factor'


        # TO EXPLORE A SINGLE SET OF VALUES FOR MODEL PARAMETERS

    return [width, height, alpha, motion, speciation, \
            seedCom, m, r, nNi, nP, nC, rmax, gmax, maintmax, dmax, \
            envgrads, rates, pmax, mmax]
