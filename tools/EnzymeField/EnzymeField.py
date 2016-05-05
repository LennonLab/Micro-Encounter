from __future__ import division
from random import sample
import numpy as np



def EnzymeField(Field, IndX, IndY, ADList, Qs, width):

    """ A function to set the initial enzyme density field

    enzymeX: X coordinates of enzyme cloud centroids
    enzymeY: Y coordinates of enzyme cloud centroids
    IndX: Individual X coordinates
    IndY: Individual Y coordinates
    ADList: A list of the metabolic states (a or d) for each individual
    Qs: A list of cell quotas for each individual

    Each individual will potentially have a local enzyme cloud

    How it works:
    1) Active individuals increase the enzyme cloud around them
    2) Enzyme clouds decay exponentially in the absence of activity

    """

    # 1.) Decay each value in the field, 0's just stay 0's

    for i, val in enumerate(Field):
        Field[i] = val/2.0  # simple half-life decay

    # 2.) Increase values in the enzyme field according to individual activity

    index = 0
    for i, Q in enumerate(Qs): # for each individual in the community

        X = int(round(IndX[i])) # get the x coordinate
        Y = int(round(IndY[i])) # get the y coordinate

        # Find the individual's location in the density field
        index = int(round(X + (Y * width)))

        # make sure an individual's location isn't out of bounds
        if index > len(Field) - 1:
            index = len(Field) - 1
        elif index < 0:
            index = 0

        # Change the field value according to individual activity and cell quota
        if ADList[i] == 'a':
            # Individual influence on the local enzyme value is a simple
            # linear relation to mean cell quota
            Field[index] += np.mean(Q)


    return Field
