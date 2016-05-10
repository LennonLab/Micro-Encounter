from __future__ import division
from random import randint
import numpy as np
import sys


def EnzymeField(Field, IndX, IndY, ADList, Qs, width):

    print 'enzyme field'

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




def consume(field, RList, R_Vals, R_IDs, R_ID, RXs, RYs, Sp_IDs, Qs, I_IDs, I_ID,
    IXs, IYs, w, h, GD, RD, DispD, GrowthList, MaintList, DispList, ADList, TLList, EVList, TrophicComplexityLevel, SpatialComplexityLevel, ResourceComplexityLevel, BiologicalComplexityLevel, enz_action = True):


    listlen = [len(Sp_IDs), len(Qs), len(I_IDs), len(IXs), len(IYs), len(GrowthList), len(MaintList), len(DispList), len(ADList), len(EVList), len(TLList)]
    if min(listlen) != max(listlen):
        print 'In consume (top)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()

    if Sp_IDs == [] or R_IDs == []:
        return [RList, R_Vals, R_IDs, R_ID, RXs, RYs, Sp_IDs, Qs]

    n = len(I_IDs)
    for ii in range(n):
        i = randint(0, n-1)

        state = ADList[i]
        if state == 'd':
            continue

        # Trophic level
        tl = TLList[i]

        x1 = IXs[i]
        y1 = IYs[i]

        Try = min([10, len(R_IDs)])
        ct = 0

        while ct < Try and R_IDs != []:
            ct += 1

            r = len(R_IDs)
            Try = min([10, r])
            j = randint(0, r-1)

            # The food
            R = RList[j]
            Rtype = R[0]

            if ResourceComplexityLevel == 1:
                tl = str(Rtype)

            if tl == Rtype: # individual is capable of consuming the resource type

                x2 = RXs[j]
                y2 = RYs[j]
                dist = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

                ind_radius = np.mean(Qs[i])
                R_val = R_Vals[j]
                res_radius = R_val

                if dist <= ind_radius + res_radius:

                    if '-' in R:
                        # useful stackoverflow answers:
                        # http://stackoverflow.com/questions/32212998/split-string-by-n-when-n-is-random
                        # http://stackoverflow.com/questions/19954593/python-checking-a-strings-first-and-last-character
                        # http://stackoverflow.com/questions/5188792/how-to-check-a-string-for-specific-characters

                        pos = randint(0, len(R))
                        if R[pos] == '-':
                            pieces = randSplit(R, randint(0,len(R)))

                            for piece in pieces:
                                switch = 'off'
                                if piece.startswith('-'):
                                    piece = piece[1:]

                                if piece.endswith('-'):
                                    piece = piece[:-1]

                                if '-' not in piece:
                                    RList[j] = piece
                                    R_Val = piece.count(Rtype)
                                    R_Vals[j] = R_Val

                                elif '-' in piece and switch == 'off':
                                    switch = 'on'
                                    RList[j] = piece
                                    val = piece.count(Rtype)
                                    R_Vals[j] = val


                                if '-' in piece and switch == 'on':
                                    RList.append(piece)
                                    R_Vals.append(val)
                                    R_IDs.append(R_ID)
                                    R_ID += 1
                                    RXs.append(x2)
                                    RYs.append(y2)
