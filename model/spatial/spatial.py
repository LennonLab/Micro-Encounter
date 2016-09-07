from __future__ import division
#import sys
from random import choice, randint
import numpy as np
import math


def distance(p0, p1):

    """ take two (x, y) tuples as parameters
    http://stackoverflow.com/questions/5407969/distance-formula-between-two-points-in-a-list"""

    return math.sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2)



def nearest_neighbor(xlist1, xlist2, ylist1, ylist2, q=1):

    """ xlist1 and ylist1 are assumed to be the lists of individual organisms
    xlist2 and ylist2 are assumed to be the lists of whatever the individual
    organisms are being measured with respect to their distance to """

    n = len(xlist1)
    r = len(xlist2)
    r = min([10, r])

    refpoints = min([10, n])
    DistList = []

    for ref in range(refpoints):

        i = randint(0, n-1)
        x1 = xlist1[i]
        y1 = ylist1[i]
        MinDist = 10000

        for j in range(r):

            x2 = xlist2[j]
            y2 = ylist2[j]
            dist = distance((x1, y1), (x2, y2))

            if dist < MinDist:
                MinDist = dist

        DistList.append(MinDist)

    return np.mean(DistList)



def morisitas(Xs, Ys, w, h):

    '''
    From http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/dispindmorisita.html

    The Morisita index of dispersion is defined as (Morisita 1959, 1962):

    Imor = n * (sum(xi^2) - sum(xi)) / (sum(xi)^2 - sum(xi))

    where xi is the count of individuals in sample i, and n is the number of samples (i = 1, 2, ..., n).
    Imor has values from 0 to n. In uniform (hyperdispersed) patterns its value falls between 0 and 1, in clumped patterns it falls between 1 and n.
    For increasing sample sizes (i.e. joining neighbouring quadrats), Imor goes to n as the quadrat size approaches clump size.
    For random patterns, Imor = 1 and counts in the samples follow Poisson frequency distribution.
    '''

    Imor = float()

    Boxes = [list([]) for _ in xrange(w*h)]

    index = 0
    for i, val in enumerate(Xs):
        X = int(round(Xs[i]))
        Y = int(round(Ys[i]))

        index = int(round(X + (Y * w)))

        if index > len(Boxes) - 1:
            index = len(Boxes) - 1
        elif index < 0:
            index = 0

        Boxes[index].append(val)

    a = 0
    b = 0

    for i, box in enumerate(Boxes):
        a += len(box)**2
        b += len(box)

    if b <= 1:
        return float('nan')
        
    Imor = len(Boxes) * (a - b)/(b**2 - b)

    return Imor





def avg_dist(xlist1, xlist2, ylist1, ylist2, q=1):

    """ xlist1 and ylist1 are assumed to be the lists of individual organisms
    xlist2 and ylist2 are assumed to be the lists of whatever the individual
    organisms are being measured with respect to their distance to """

    nmax = len(xlist1)
    rmax = len(xlist2)

    refpoints = min([100, nmax])
    dist = []

    for n in range(refpoints):
        for j in range(q):

            i1 = choice(range(nmax))
            x1 = xlist1[i1]
            y1 = ylist1[i1]

            i2 = choice(range(rmax))
            x2 = xlist2[i2]
            y2 = ylist2[i2]

            dist.append(distance((x1, y1), (x2, y2)))

    return np.mean(dist)


"""
def nearest_neighbor_dist(xlist1, xlist2, ylist1, ylist2, q=1):

    nmax = len(xlist1)
    refpoints = min([100, nmax])

    for n in range(refpoints):
        for j in range(q):

            i = choice(range(nmax))
            x1 = xlist1[i]
            x2 = xlist2[i]
            y1 = ylist1[i]
            y2 = ylist2[i]


    return
"""
