# -*- coding: utf-8 -*-
from __future__ import division
from random import randint, choice
import numpy as np
import math
import sys
import os


def check_list_lengths(lists, function, ct):
    lengths = []
    for lst in lists:
        lengths.append(len(lst))

    if min(lengths) != max(lengths):
        print ct, function, 'min(listlen) != max(listlen)'
        print lengths
        sys.exit()


def GetRAD(vector):
    RAD = []
    unique = list(set(vector))
    for val in unique:
        RAD.append(vector.count(val)) # the abundance of each Sp_

    return RAD, unique # the rad and the specieslist


def per_capita(ValDict, SpeciesIDs):
    species, vals = [ValDict.keys(), ValDict.values()]
    for i, sp in enumerate(species):
        sp_val = SpeciesIDs.count(sp) * vals[i]

    return sp_val/len(SpeciesIDs) # the rad and the specieslist



def ResIn(ResLists, RID, params, ct, ComplexityLevels):

    Rvals, RX, RY, RIDs = ResLists
    width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std = params
    SC, TC, RC = ComplexityLevels

    Ymean = float(np.random.uniform(0.001*height, 0.999*height))
    Xmean = float(np.random.uniform(0.001*width, 0.999*width))

    for i in range(r):
        RIDs.append(RID)
        RID += 1

        if RC == 1:
            #http://stackoverflow.com/questions/1957273/how-do-i-generate-a-random-string-of-length-x-a-z-only-in-python
            # all resource particles are of the same type
            ri = choice(['a', 'aa', 'aaaa'])
            if ri == 'a': n = randint(1, 40)
            elif ri == 'aa': n = randint(1, 20)
            elif ri == 'aaaa': n = randint(1, 10)
            ri = ''.join(ri for _ in xrange(n))

        elif RC == 2:
            # Three types of resource particles.
            ri = choice(['a', 'bb', 'cccc'])
            if ri == 'a': n = randint(1, 40)
            elif ri == 'bb': n = randint(1, 20)
            elif ri == 'cccc': n = randint(1, 10)

            ri = ''.join(ri for _ in xrange(n))

        elif RC == 3:
            # Three types of resource particles that vary in the effort/energy needed to break them down
            ri = choice(['a-', 'bb-', 'c-cc-c-'])
            if ri == 'a-': n = randint(1, 40)
            elif ri == 'bb-': n = randint(1, 20)
            elif ri == 'c-cc-c-': n = randint(1, 10)

            ri = ''.join(ri for _ in xrange(n))

        Rvals.append(ri)

        if SC == 1:
            RY.append(float(np.random.uniform(0.001*height, 0.999*height)))
            RX.append(float(np.random.uniform(0.001*width, 0.999*width)))

        elif SC <= 3:
            vals = np.random.normal([Ymean, Xmean], std)

            if vals[0] < 0.001*height: vals[0] = 0.001*height
            elif vals[0] > 0.999*height: vals[0] = 0.999*height

            if vals[1] < 0.001*width: vals[1] = 0.001*width
            elif vals[1] > 0.999*width: vals[1] = 0.999*width

            RY.append(vals[0])
            RX.append(vals[1])

    ResLists = Rvals, RX, RY, RIDs

    return [ResLists, RID]



def immigration(IndLists, IndDicts, IndID, params, ct, ComplexityLevels):

    SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList = IndLists
    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = IndDicts
    width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std = params
    SC, TC, RC = ComplexityLevels

    Ymean = float(np.random.uniform(0.001*height, 0.999*height))
    Xmean = float(np.random.uniform(0.001*width, 0.999*width))

    check_list_lengths(IndLists, 'immigration', ct)

    if ct > 1: seedCom = 1
    if ct == 1: SC = 1

    for m in range(seedCom):
        prop = np.random.randint(1, 100)
        SpeciesIDs.append(prop)

        if SC <= 3:
            IndY.append(float(np.random.uniform(0.001*height, 0.999*height)))
            IndX.append(float(np.random.uniform(0.001*width, 0.999*width)))

        elif SC <= 3:
            vals = np.random.normal([Ymean, Xmean], std)

            if vals[0] < 0.001*height: vals[0] = 0.001*height
            elif vals[0] > 0.999*height: vals[0] = 0.999*height

            if vals[1] < 0.001*width: vals[1] = 0.001*width
            elif vals[1] > 0.999*width: vals[1] = 0.999*width

            IndY.append(vals[0])
            IndX.append(vals[1])

        IndIDs.append(IndID)
        IndID += 1

        Q = float(np.random.uniform(0.1, 0.1))
        CellQuotas.append(Q)

        if prop not in GrowthDict:

            # species trophic level
            TrophicDict[prop] = choice(['a', 'b', 'c'])

            # species growth rate
            GrowthDict[prop] = np.random.uniform(gmax/10, gmax)

            # species maintenance
            MaintDict[prop] = np.random.uniform(maintmax/10, maintmax)

            # species maintenance reduction factor
            MainFactorDict[prop] = np.random.uniform(10, mmax)

            # species resuscitation factor
            RPFDict[prop] = np.random.uniform(pmax/10, pmax)

            # species active dispersal rate
            DispDict[prop] = np.random.uniform(dmax/10, dmax)

        ADList.append('a')

    IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
    IndDicts = GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict

    return [IndLists, IndDicts, IndID]




def maintenance(IndLists, IndDicts, ResLists, RID, ComplexityLevels, numDead):

    SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList = IndLists
    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = IndDicts
    Rvals, RX, RY, RIDs = ResLists
    SC, TC, RC = ComplexityLevels

    n = len(IndIDs)
    for j in range(n):

        i = randint(0, len(IndIDs)-1)
        Q = CellQuotas[i]
        spID = SpeciesIDs[i]
        mfd = MainFactorDict[spID]  # multiplicative maintenance factor (is greater than 1)
        maint = MaintDict[spID]
        state = ADList[i]

        if state == 'd': maint = maint/mfd

        Q -= maint
        if Q <= maint:   # starved

            if TC >= 3:
                r = 'd'
                n = 1
                if RC == 3:
                    r, n = choice([['d-',4], ['dd-',2], ['d-dd-d-',4]])

                r = ''.join(r for _ in xrange(n))
                Rvals.append(r)
                RX.append(IndX[i])
                RY.append(IndY[i])
                RIDs.append(RID)
                RID += 1

            CellQuotas.pop(i)
            SpeciesIDs.pop(i)
            IndIDs.pop(i)
            IndX.pop(i)
            IndY.pop(i)
            ADList.pop(i)
            numDead += 1

        else: CellQuotas[i] = Q

    IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
    IndDicts = GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict

    return [IndLists, ResLists, RID, numDead]





def transition(IndLists, IndDicts, ResLists, RID, ComplexityLevels, numDead):

    SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList, = IndLists
    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = IndDicts
    Rvals, RX, RY, RIDs = ResLists
    SC, TC, RC = ComplexityLevels

    n = len(IndIDs)
    for j in range(n):

        i = randint(0, len(IndIDs)-1)
        spID = SpeciesIDs[i]
        state = ADList[i]
        Q = CellQuotas[i]

        mfd = MainFactorDict[spID]  # multiplicative maintenance factor (is greater than 1)
        maint = MaintDict[spID]

        if state == 'd':
            x = np.random.binomial(1, RPFDict[spID]) # make this probability a randomly chosen variable
            if x == 1:
                ADList[i] = 'a'
                Q -= maint/mfd

        if state == 'a':
            if Q <= maint:  # go dormant
                ADList[i] = 'd'

            if Q < 0.0:

                if TC >= 3:
                    r = 'd'
                    n = 1

                    if RC == 3:
                        r, n = choice([['d-',4], ['dd-',2], ['d-dd-d-',4]])

                    r = ''.join(r for _ in xrange(n))
                    Rvals.append(r)

                    RX.append(RX[i])
                    RY.append(RY[i])
                    RIDs.append(RID)
                    RID += 1

                CellQuotas.pop(i)
                SpeciesIDs.pop(i)
                IndIDs.pop(i)
                IndX.pop(i)
                IndY.pop(i)
                ADList.pop(i)
                numDead += 1
                continue


    IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
    ResLists = Rvals, RX, RY, RIDs

    return [IndLists, ResLists, RID, numDead]




def consume(IndLists, IndDicts, IndID, ResLists, RID, params, ComplexityLevels):

    SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList, = IndLists
    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = IndDicts
    width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std = params
    Rvals, RX, RY, RIDs = ResLists
    SC, TC, RC = ComplexityLevels

    encounters = 0
    n = len(IndIDs)
    for ii in range(n):
        i = randint(0, n-1)

        # Trophic level
        spID = SpeciesIDs[i]
        tl = TrophicDict[spID]
        Q = CellQuotas[i]

        mfd = MainFactorDict[spID]
        x1 = IndX[i]
        y1 = IndY[i]

        Try = min([40, len(RIDs)])
        ct = 0

        while ct < Try and RIDs != []:
            ct += 1

            r = len(RIDs)
            Try = min([40, r])
            j = randint(0, r-1)
            R = Rvals[j]

            if R.startswith('-'): R = R[1:]
            if R.endswith('-'): R = R[:-1]

            Rtype = R[0]
            if Rtype == '-':
                print 'Rtype is a hyphen'
                sys.exit()

            if Rtype == 'd': tl = str(Rtype)

            if tl == Rtype: # is capable of consuming the resource type
                x2 = RX[j]
                y2 = RY[j]
                dist = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

                ind_radius = np.mean(CellQuotas[i])
                Rval = R.count(Rtype)
                res_radius = Rval/20.0
                tl = TrophicDict[spID] # need to reset this in case of no match

                if dist <= ind_radius + res_radius:
                    if ADList[i] == 'd':
                        ADList[i] = 'a'
                        MaintDict[spID] = MaintDict[spID]*mfd # metabolic maintenance increases

                    eat = 'no'
                    if Rtype == 'd':
                        eat = 'yes'

                    ct = Try
                    if '-' in R:
                        pos = randint(1, len(R)-1)

                        if R[pos] == '-':
                            pieces = [R[:pos], R[pos:]]

                            for index, piece in enumerate(pieces):
                                if piece.startswith('-'): piece = piece[1:]

                                if piece.endswith('-'): piece = piece[:-1]

                            # Split the molecule but do not consume
                            Rvals[j] = pieces[0]
                            Rvals.append(pieces[1])
                            RIDs.append(RID)
                            RID += 1
                            RX.append(x2)
                            RY.append(y2)


                    elif '-' not in R:
                        eat = 'yes'

                        if len(R) > 1:
                            pieces = [R[:1], R[1:]]
                            Rvals[j] = pieces[0]
                            Rvals.append(pieces[1])
                            RIDs.append(RID)
                            RID += 1
                            RX.append(x2)
                            RY.append(y2)

                    if eat == 'yes':
                        encounters += 1
                        mu = GrowthDict[spID]
                        # Increase cell quota & decrease resource particle size
                        Q += (mu * Q)
                        if Q > 1.0: Q = 1.0

                        levels = [2, 4]
                        if TC in levels and Rtype != 'c':
                            # b's are by-products of consuming a's
                            # c's are by-products of consuming b's

                            if Rtype == 'a': R = R.replace('a', 'b')
                            if Rtype == 'b': R = R.replace('b', 'c')

                            Rvals[j] = R
                            RIDs[j] = RID
                            RID += 1

                        else:
                            # no resources produced as by-products
                            Rvals.pop(j)
                            RIDs.pop(j)
                            RX.pop(j)
                            RY.pop(j)


                        CellQuotas[i] = Q

    IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
    ResLists = Rvals, RX, RY, RIDs

    return [ResLists, RID, IndLists, encounters]





def reproduce(IndLists, IndDicts, IndID, ResLists, RID, params, ComplexityLevels, numDead):

    SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList, = IndLists
    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = IndDicts
    width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std = params
    Rvals, RX, RY, RIDs = ResLists
    SC, TC, RC = ComplexityLevels

    n = len(IndIDs)
    for j in range(n):

        i = randint(0, len(IndIDs)-1)
        state = ADList[i]
        if state == 'd': continue

        Q = CellQuotas[i]
        pq = float(np.mean(Q))

        if Q < 0.0:

            if TC >= 3:
                r = 'd'
                n = 1

                if RC == 3:
                    r, n = choice([['d-',4], ['dd-',2], ['d-dd-d-',4]])

                r = ''.join(r for _ in xrange(n))
                Rvals.append(r)
                RX.append(RX[i])
                RY.append(RY[i])
                RIDs.append(RID)
                RID += 1

            CellQuotas.pop(i)
            SpeciesIDs.pop(i)
            IndIDs.pop(i)
            IndX.pop(i)
            IndY.pop(i)
            ADList.pop(i)
            numDead += 1
            continue

        p = np.random.binomial(1, pq)
        if p == 1: # individual is large enough to reproduce

            spID = SpeciesIDs[i]
            X = IndX[i]
            Y = IndY[i]

            p = 1
            if p == 1: # the environment is suitable for reproduction

                CellQuotas[i] = Q/2
                CellQuotas.append(Q/2)

                IndID += 1
                IndIDs.append(IndID)

                SpeciesIDs.append(spID)
                ADList.append('a')

                newX = float(np.random.uniform(X-0.1, X+0.1, 1))
                if 0.1 > newX: newX = 0
                if newX > width - 0.1: newX = width - 0.1
                IndX.append(newX)

                newY = float(np.random.uniform(Y-0.1, Y+0.1, 1))
                if 0.1 > newY: newY = 0
                elif newY > height: newY = height - 0.1
                IndY.append(newY)


    IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
    ResLists = Rvals, RX, RY, RIDs

    return [IndLists, IndDicts, IndID, ResLists, RID, numDead]




def res_dispersal(ResLists, RID, params, ct, ComplexityLevels):

    Rvals, RX, RY, RIDs = ResLists
    width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std = params

    n = len(RIDs)
    for j in range(n):

        i = randint(0, len(RIDs)-1)

        RY[i] = float(np.random.uniform(0.001*height, 0.999*height))
        RX[i] = float(np.random.uniform(0.001*width, 0.999*width))

    return [ResLists, RID]



def dispersal(IndLists, IndDicts, IndID, ResLists, RID, params, ComplexityLevels, numDead):

    SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList = IndLists
    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = IndDicts
    width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std = params
    Rvals, RX, RY, RIDs = ResLists
    SC, TC, RC = ComplexityLevels

    n = len(IndIDs)
    for j in range(n):

        i = randint(0, len(IndIDs)-1)
        spID = SpeciesIDs[i]
        X = IndX[i]
        Y = IndY[i]


        if SC == 1:
            IndY[i] = float(np.random.uniform(0.001*height, 0.999*height))
            IndX[i] = float(np.random.uniform(0.001*width, 0.999*width))

        state = ADList[i]
        if state == 'd': continue

        if SC == 2:
            dist = DispDict[spID]

            Q = CellQuotas[i]
            if Q >= MaintDict[spID]*dist:

                # A cost for active dispersal
                Q -= MaintDict[spID]*dist

                if Q < 0.0:

                    if TC >= 3:
                        r = 'd'
                        n = 1

                        if RC == 3:
                            r, n = choice([['d-',4], ['dd-',2], ['d-dd-d-',4]])

                        r = ''.join(r for _ in xrange(n))
                        Rvals.append(r)

                        RX.append(RX[i])
                        RY.append(RY[i])
                        RIDs.append(RID)
                        RID += 1

                    CellQuotas.pop(i)
                    SpeciesIDs.pop(i)
                    IndIDs.pop(i)
                    IndX.pop(i)
                    IndY.pop(i)
                    ADList.pop(i)
                    numDead += 1
                    continue

                CellQuotas[i] = Q
                vd = choice([-1, 1])
                hd = choice([-1, 1])

                X += hd*dist
                Y += vd*dist

                if X > 0.999 * width: X = width * 0.999
                elif X < 0.001 * width: X = 0.001 * width
                if Y > 0.999 * height: Y = height * 0.999
                elif Y < 0.001 * height: Y = 0.001 * height

                IndX[i] = X
                IndY[i] = Y

    IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
    ResLists = Rvals, RX, RY, RIDs

    return [IndLists, IndDicts, IndID, ResLists, RID, numDead]




def chemotaxis(IndLists, IndDicts, IndID, ResLists, RID, params, ComplexityLevels, numDead):

    SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList = IndLists
    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = IndDicts
    width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std = params
    Rvals, RX, RY, RIDs = ResLists
    SC, TC, RC = ComplexityLevels

    n = len(IndIDs)
    for ii in range(n):

        n = len(IndIDs)
        i = randint(0, n-1)
        state = ADList[i]

        if state == 'd': continue

        spID = SpeciesIDs[i]
        disp = DispDict[spID]
        tl = TrophicDict[spID] # Trophic level
        x1 = IndX[i]
        y1 = IndY[i]

        Try = min([40, len(RIDs)])
        ct = 0

        minDist = 100
        targetX = 0
        targetY = 0
        dist = 0

        while ct < Try and RIDs is not []:

            ct += 1
            r = len(RIDs)
            Try = min([40, r])
            j = randint(0, r-1)

            # The food
            R = Rvals[j]
            if R.startswith('-'): R = R[1:]
            if R.endswith('-'): R = R[:-1]
            Rtype = R[0]


            if tl == Rtype: # is capable of consuming the resource type
                x2 = RX[j]
                y2 = RY[j]
                dist = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

                if dist < minDist:
                    minDist = dist
                    targetX = RX[j]
                    targetY = RY[j]


        # A cost for active dispersal
        Q = CellQuotas[i]
        maint = MaintDict[spID]
        if Q - maint*disp > maint:
            Q -= maint*disp
        else:
            dist = (Q - maint)/maint
            Q -= maint*disp

        if Q < 0.0:

            if TC >= 3:
                r = 'd'
                n = 1

                if RC == 3:
                    r, n = choice([['d-',4], ['dd-',2], ['d-dd-d-',4]])

                r = ''.join(r for _ in xrange(n))
                Rvals.append(r)
                RX.append(IndX[i])
                RY.append(IndY[i])
                RIDs.append(RID)
                RID += 1

            CellQuotas.pop(i)
            SpeciesIDs.pop(i)
            IndIDs.pop(i)
            IndX.pop(i)
            IndY.pop(i)
            ADList.pop(i)
            numDead += 1
            continue

        CellQuotas[i] = Q

        if x1 > targetX: x1 -= dist*disp
        elif x1 < targetX: x1 += dist*disp
        if y1 > targetY: y1 -= dist*disp
        elif y1 < targetY: y1 += dist*disp

        if x1 > 0.999 * width: x1 = width * 0.999
        elif x1 < 0.001 * width: x1 = 0.001 * width
        if y1 > 0.999 * height: y1 = height * 0.999
        elif y1 < 0.001 * height: y1 = 0.001 * height

        IndX[i] = x1
        IndY[i] = y1

    IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
    ResLists = Rvals, RX, RY, RIDs

    return [ResLists, RID, IndLists, IndDicts, IndID, numDead]
