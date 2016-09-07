# -*- coding: utf-8 -*-
from __future__ import division
from random import randint, choice
import numpy as np
import math
import string
import sys
import os

mydir = os.path.expanduser("~/GitHub/Micro-Encounter")
sys.path.append(mydir + "/model/metrics")
import metrics


def checkQ(Q, line):
    if Q < 0:
        print 'line',line,': error: Q < 0:', Q
        sys.exit()
    return


def dead(IndLists, ResLists, RID, i, Q, TC, RC):

    SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList = IndLists
    Rvals, Rtypes, RX, RY, RIDs = ResLists

    if 'scavenging' in TC and 'none' not in TC:
        Rvals.append(Q)
        Rtypes.append('dead')
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

    IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
    ResLists = Rvals, Rtypes, RX, RY, RIDs
    return [IndLists, ResLists, RID]



def ResIn(ResLists, ResDict, RID, params, ct, ComplexityLevels):

    Rvals, Rtypes, RX, RY, RIDs = ResLists
    width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std = params
    SC, TC, RC = ComplexityLevels

    Ymean = np.random.uniform(0.001*height, 0.999*height)
    Xmean = np.random.uniform(0.001*width, 0.999*width)
    sizemax = 100

    for i in range(r):
        RIDs.append(RID)
        RID += 1
        Rvals.append(randint(1, sizemax))
        pmin, pmax = float(), float() # default is 'simple'
        ri = str()

        if '-polyculture-' in RC:
            ri = choice(['a','b','c','d','e'])
            #ri = choice(string.lowercase) # default is polyculture
        elif '-monoculture-' in RC:
            ri = 'a'

        if '-simple-' in RC:
            pmin, pmax = 1.0, 1.0
        elif '-lockandkey-' in RC:
            pmin, pmax = 0.05, 0.1

        if ri not in ResDict:
            ResDict[ri] = round(np.random.uniform(pmin, pmax), 3)

        Rtypes.append(ri)

        if '-aggregated-' in SC:
            vals = np.random.normal([Ymean, Xmean], std)
            RY.append(vals[0])
            RX.append(vals[1])

        elif '-random-' in SC:
            RY.append(np.random.uniform(0.001*height, 0.999*height))
            RX.append(np.random.uniform(0.001*width, 0.999*width))

    ResLists = Rvals, Rtypes, RX, RY, RIDs

    return [ResLists, ResDict, RID]



def immigration(IndLists, IndDicts, IndID, params, ct, ComplexityLevels):

    SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList = IndLists
    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = IndDicts
    width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std = params
    SC, TC, RC = ComplexityLevels

    Ymean = np.random.uniform(0.001*height, 0.999*height)
    Xmean = np.random.uniform(0.001*width, 0.999*width)

    if ct > 1: seedCom = 1
    for m in range(seedCom):
        prop = np.random.randint(1, 100)
        SpeciesIDs.append(prop)

        CellQuotas.append(np.random.uniform(100, 100))
        if '-aggregated-' in SC:
            vals = np.random.normal([Ymean, Xmean], std)
            IndY.append(vals[0])
            IndX.append(vals[1])

        elif '-random-' in SC:
            IndY.append(np.random.uniform(0.001*height, 0.999*height))
            IndX.append(np.random.uniform(0.001*width, 0.999*width))

        IndIDs.append(IndID)
        IndID += 1

        pred = np.random.binomial(1, 0.01)

        if prop not in GrowthDict:
            TrophicDict[prop] = {'resource': str(), 'crossfeed_to': str(), 'prey_of': int(), 'predator_of': int(), 'mutualist_with': int()}

            if '-monoculture-' in RC:
                TrophicDict[prop]['resource'] = 'a'

            elif '-polyculture-' in RC:
                TrophicDict[prop]['resource'] = choice(['a','b','c','d','e'])
                #TrophicDict[prop]['resource'] = choice(string.lowercase)

            if '-crossfeeding-' in TC:
                TrophicDict[prop]['crossfeeds_to'] = choice(['a','b','c','d','e'])
                #TrophicDict[prop]['crossfeeds_to'] = choice(string.lowercase)

            if '-predprey-' in TC and pred == 1:
                prey_of = int()
                pred_of = int()

                while prey_of is not prop and pred_of is not prop:
                    prey_of = np.random.randint(1, 100)
                    pred_of = np.random.randint(1, 100)
                    if prey_of != pred_of and prop not in [prey_of, pred_of]:
                        break

                TrophicDict[prop]['prey_of'] = prey_of
                TrophicDict[prop]['predator_of'] = pred_of

            # species growth rate
            GrowthDict[prop] = np.random.uniform(gmax/1, gmax)

            # species maintenance
            MaintDict[prop] = randint(1, maintmax)

            # species maintenance reduction factor
            MainFactorDict[prop] = randint(1, mmax)

            # species resuscitation factor
            RPFDict[prop] = np.random.uniform(pmax/1, pmax)

            # species active dispersal rate
            DispDict[prop] = np.random.uniform(dmax/1, dmax)

        ADList.append('a')

    metrics.check_list_lengths(IndLists, 'immigration', ct)

    IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
    IndDicts = GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict
    return [IndLists, IndDicts, IndID]




def maintenance(IndLists, IndDicts, ResLists, ResDict, RID, ComplexityLevels, numDead):

    SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList = IndLists
    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = IndDicts
    Rvals, Rtypes, RX, RY, RIDs = ResLists
    SC, TC, RC = ComplexityLevels

    n = len(IndIDs)
    if n == 0: return [IndLists, ResLists, ResDict, RID, numDead]

    for j in range(n):
        i = randint(0, len(IndIDs)-1)
        Q = CellQuotas[i]
        spID = SpeciesIDs[i]
        mfd = MainFactorDict[spID]
        maint = MaintDict[spID]
        state = ADList[i]

        if state == 'd':
            maint = maint/mfd

        if Q >= maint/mfd:
            #checkQ(Q, 182)
            Q -= maint/mfd
            #checkQ(Q, 184)
            CellQuotas[i] = Q

        else:
            IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
            ResLists = Rvals, Rtypes, RX, RY, RIDs
            IndLists, ResLists, RID = dead(IndLists, ResLists, RID, i, Q, TC, RC)
            numDead += 1

    IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
    IndDicts = GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict
    return [IndLists, ResLists, ResDict, RID, numDead]




def transition(IndLists, IndDicts, ResLists, ResDict, RID, ComplexityLevels, numDead):

    SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList = IndLists
    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = IndDicts
    Rvals, Rtypes, RX, RY, RIDs = ResLists
    SC, TC, RC = ComplexityLevels

    n = len(IndIDs)
    if n == 0:
        return [IndLists, ResLists, ResDict, RID, numDead]

    for j in range(n):
        i = randint(0, len(IndIDs)-1)
        spID = SpeciesIDs[i]
        state = ADList[i]
        Q = CellQuotas[i]

        mfd = MainFactorDict[spID]
        maint = MaintDict[spID]

        if Q < maint/mfd:
            IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
            ResLists = Rvals, Rtypes, RX, RY, RIDs
            IndLists, ResLists, RID = dead(IndLists, ResLists, RID, i, Q, TC, RC)
            numDead += 1
            continue

        if state == 'd':
            x = np.random.binomial(1, RPFDict[spID])
            if x == 1: ADList[i] = 'a'

        if state == 'a':
            if Q <= maint: ADList[i] = 'd'

    IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
    ResLists = Rvals, Rtypes, RX, RY, RIDs
    return [IndLists, ResLists, ResDict, RID, numDead]




def consume(IndLists, IndDicts, IndID, ResLists, ResDict, RID, params, ComplexityLevels, numDead):

    SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList = IndLists
    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = IndDicts
    width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std = params
    Rvals, Rtypes, RX, RY, RIDs = ResLists
    SC, TC, RC = ComplexityLevels

    maxQ = 1000
    encounters = 0
    n = len(IndIDs)
    r = len(RIDs)

    if n == 0 and r == 0:
        return [ResLists, ResDict, RID, IndLists, encounters, numDead]

    elif r == 0 and n > 0:
        if '-predprey-' in TC and 'none' not in TC:
            ResLists, ResDict, RID, IndLists, encounters, numDead = predation(IndLists, IndDicts, IndID, ResLists, ResDict, RID, params, ComplexityLevels, encounters, numDead)
            SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList = IndLists
            Rvals, Rtypes, RX, RY, RIDs = ResLists

        return [ResLists, ResDict, RID, IndLists, encounters, numDead]

    n = len(IndIDs)
    r = len(RIDs)

    for ii in range(n):
        i = randint(0, n-1)

        # Trophic level
        spID = SpeciesIDs[i]
        Q = CellQuotas[i]
        x1 = IndX[i]
        y1 = IndY[i]

        restype = TrophicDict[spID]['resource']
        res_indices = [xi for xi, x in enumerate(Rtypes) if x == restype or x == 'dead']

        if len(res_indices) == 0: continue

        Try = min([20, len(res_indices)])
        ct = 0

        while ct < Try and RIDs != []:
            ct += 1

            r = len(res_indices)
            Try = min([20, r])
            j = randint(0, r-1)
            Rtype = Rtypes[j]

            if restype == Rtype or Rtype == 'dead': # is capable of consuming the resource type
                x2 = RX[j]
                y2 = RY[j]
                dist = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

                Rval = Rvals[j]
                if dist <= Q + Rval:

                    if ADList[i] == 'd':
                        ADList[i] = 'a'

                    ct = Try
                    p = ResDict[Rtype]
                    k = np.random.binomial(1, p)
                    #print k, p

                    if k == 1:
                        encounters += 1

                        p2 = np.random.uniform(0, 1.0)
                        res_piece1 = Rval * p2
                        res_piece2 = Rval - res_piece1

                        Rvals[j] = float(res_piece1)
                        Rvals.append(res_piece2)
                        RIDs.append(RID)
                        Rtypes.append(Rtype)
                        RID += 1
                        RX.append(x2)
                        RY.append(y2)

                        mu = GrowthDict[spID] * Q
                        Rval = Rvals[j]

                        if mu >= Rval:
                            Q += Rval
                            Rval = 0.0

                        elif mu < Rval:
                            Q += mu
                            Rval -= mu

                        if Q > maxQ:
                            Rval += maxQ - Q
                            Q = maxQ

                        if Rval == 0.0:
                            Rvals.pop(j)
                            RIDs.pop(j)
                            Rtypes.pop(j)
                            RX.pop(j)
                            RY.pop(j)

                        else:
                            if 'polyculture' in RC:
                                if 'crossfeeding' in TC and 'none' not in TC:
                                    Rtype = TrophicDict[spID]['crossfeeds_to']

                                    if 'simple' in RC:
                                        pmin, pmax = 1.0, 1.0
                                    if 'lockandkey' in RC:
                                        pmin, pmax = 0.05, 0.1

                                    if Rtype not in ResDict:
                                        ResDict[Rtype] = round(np.random.uniform(pmin, pmax), 3)

                            Rvals[j] = Rval
                            Rtypes[j] = Rtype
                            RIDs[j] = RID
                            RID += 1

                        CellQuotas[i] = Q


    if 'predprey' in TC and 'none' not in TC:
        ResLists, ResDict, RID, IndLists, encounters, numDead = predation(IndLists, IndDicts, IndID, ResLists, ResDict, RID, params, ComplexityLevels, encounters, numDead)
        SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList = IndLists
        Rvals, Rtypes, RX, RY, RIDs = ResLists

    IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
    ResLists = Rvals, Rtypes, RX, RY, RIDs
    return [ResLists, ResDict, RID, IndLists, encounters, numDead]




def predation(IndLists, IndDicts, IndID, ResLists, ResDict, RID, params, ComplexityLevels, encounters, numDead):

    SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList = IndLists
    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = IndDicts
    Rvals, Rtypes, RX, RY, RIDs = ResLists
    width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std = params
    SC, TC, RC = ComplexityLevels

    maxQ = 1000
    n = len(IndIDs)

    for ii in range(n):
        i = randint(0, n-1)

        pred = SpeciesIDs[i]
        prey = TrophicDict[pred]['predator_of']
        if prey == 0:
            continue

        predQ = CellQuotas[i]
        predx = IndX[i]
        predy = IndY[i]

        prey_indices = [i for i, x in enumerate(SpeciesIDs) if x == prey]
        Try = min([20, len(prey_indices)])

        ct = 0
        while ct < Try and IndX != []:
            ct += 1

            r = len(prey_indices)
            Try = min([20, r])
            j = randint(0, r-1)

            preyx = IndX[j]
            preyy = IndY[j]
            dist = math.sqrt((predx - preyx)**2 + (predy - preyy)**2)

            preyQ = CellQuotas[j]

            if dist <= predQ + preyQ:
                if ADList[i] == 'd':
                    ADList[i] = 'a'

                mu = GrowthDict[pred]
                if mu >= preyQ:
                    predQ += preyQ
                    preyQ = 0

                elif mu < preyQ:
                    predQ += mu
                    preyQ -= mu

                if predQ > maxQ:
                    preyQ += maxQ - predQ
                    predQ = maxQ

                IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
                ResLists = Rvals, Rtypes, RX, RY, RIDs
                IndLists, ResLists, RID = dead(IndLists, ResLists, RID, j, preyQ, TC, RC)
                numDead += 1
                continue

                CellQuotas[i] = int(predQ)

    IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
    ResLists = Rvals, Rtypes, RX, RY, RIDs
    return [ResLists, ResDict, RID, IndLists, encounters, numDead]



def reproduce(IndLists, IndDicts, IndID, ResLists, ResDict, RID, params, ComplexityLevels, numDead):

    SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList = IndLists
    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = IndDicts
    width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std = params
    Rvals, Rtypes, RX, RY, RIDs = ResLists
    SC, TC, RC = ComplexityLevels

    maxQ = 1000
    PRODI = 0
    n = len(IndIDs)
    if n == 0:
        return [PRODI, IndLists, IndDicts, IndID, ResLists, ResDict, RID, numDead]

    for j in range(n):
        i = randint(0, len(IndIDs)-1)
        state = ADList[i]
        spID = SpeciesIDs[i]
        maint = MaintDict[spID]
        mfd = MainFactorDict[spID]

        if state == 'd': continue

        Q = CellQuotas[i]
        #checkQ(Q, 463)

        if Q < maint/mfd:
            IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
            ResLists = Rvals, Rtypes, RX, RY, RIDs
            IndLists, ResLists, RID = dead(IndLists, ResLists, RID, i, Q, TC, RC)
            numDead += 1
            continue

        p = np.random.binomial(1, Q/maxQ)

        if p == 1: # individual is large enough to reproduce
            PRODI += 1
            spID = SpeciesIDs[i]
            X = IndX[i]
            Y = IndY[i]

            CellQuotas[i] = Q/2
            CellQuotas.append(Q/2)
            IndID += 1
            IndIDs.append(IndID)
            SpeciesIDs.append(spID)
            ADList.append('a')

            newX = float(np.random.uniform(X-0.1, X+0.1, 1))
            IndX.append(newX)
            newY = float(np.random.uniform(Y-0.1, Y+0.1, 1))
            IndY.append(newY)

    IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
    ResLists = Rvals, Rtypes, RX, RY, RIDs
    return [PRODI, IndLists, IndDicts, IndID, ResLists, ResDict, RID, numDead]




def res_dispersal(ResLists, ResDict, RID, params, ct, ComplexityLevels):

    Rvals, Rtypes, RX, RY, RIDs = ResLists
    width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std = params
    SC, TC, RC = ComplexityLevels

    n = len(RIDs)
    for j in range(n):
        i = randint(0, len(RIDs)-1)

        if 'wellmixed' in SC:
            RY[i] = float(np.random.uniform(0.001*height, 0.999*height))
            RX[i] = float(np.random.uniform(0.001*width, 0.999*width))
        elif 'brownian' in SC:
            RY[i] += choice([-1,1]) * 0.01
            RX[i] += choice([-1,1]) * 0.01

    return [ResLists, ResDict, RID]



def dispersal(IndLists, IndDicts, IndID, ResLists, ResDict, RID, params, ComplexityLevels, numDead):

    SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList = IndLists
    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = IndDicts
    width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std = params
    Rvals, Rtypes, RX, RY, RIDs = ResLists
    SC, TC, RC = ComplexityLevels

    n = len(IndIDs)
    if n == 0:
        return [IndLists, IndDicts, IndID, ResLists, ResDict, RID, numDead]

    if 'chemotaxis' in SC or 'pred-prey' in TC:
        fxn = 'chemotaxis'
        IndLists, IndDicts, IndID, ResLists, ResDict, RID, numDead = chemotaxis(IndLists, IndDicts, IndID, ResLists, ResDict, RID, params, ComplexityLevels, numDead, fxn)
        fxn = 'pred-prey'
        return chemotaxis(IndLists, IndDicts, IndID, ResLists, ResDict, RID, params, ComplexityLevels, numDead, fxn)

    if 'wellmixed' in SC or 'randwalk' in SC:
        for j in range(n):
            i = randint(0, len(IndIDs)-1)
            spID = SpeciesIDs[i]
            disp = DispDict[spID]

            if 'wellmixed' in SC:
                IndY[i] = float(np.random.uniform(0.001*height, 0.999*height))
                IndX[i] = float(np.random.uniform(0.001*width, 0.999*width))

            elif 'randwalk' in SC:
                IndY[i] += choice([-1,1]) * disp
                IndX[i] += choice([-1,1]) * disp

                # A cost for active dispersal
                Q = CellQuotas[i]
                maint = MaintDict[spID]
                mfd = MainFactorDict[spID]
                Q -= maint*disp

                if Q < maint/mfd:
                    IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
                    ResLists = Rvals, Rtypes, RX, RY, RIDs
                    IndLists, ResLists, RID = dead(IndLists, ResLists, RID, i, Q, TC, RC)
                    numDead += 1

                else: CellQuotas[i] = Q

    IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
    ResLists = Rvals, Rtypes, RX, RY, RIDs
    return [IndLists, IndDicts, IndID, ResLists, ResDict, RID, numDead]




def chemotaxis(IndLists, IndDicts, IndID, ResLists, ResDict, RID, params, ComplexityLevels, numDead, fxn):

    SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList = IndLists
    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = IndDicts
    width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std = params
    Rvals, Rtypes, RX, RY, RIDs = ResLists
    SC, TC, RC = ComplexityLevels

    n = len(IndIDs)
    if n == 0: return [IndLists, IndDicts, IndID, ResLists, ResDict, RID, numDead]

    for ii in range(n):
        n = len(IndIDs)
        i = randint(0, n-1)
        state = ADList[i]

        if state == 'd': continue

        spID = SpeciesIDs[i]
        disp = DispDict[spID]

        if fxn == 'chemotaxis':
            restype = TrophicDict[spID]['resource']
        elif fxn == 'pred-prey':
            restype = TrophicDict[spID]['predator_of']

        x1 = IndX[i]
        y1 = IndY[i]

        res_indices = []
        if fxn == 'chemotaxis':
            res_indices = [xi for xi, x in enumerate(Rtypes) if x == restype]
        elif fxn == 'pred-prey':
            res_indices = [xi for xi, x in enumerate(SpeciesIDs) if x == restype or x == 'dead']

        if len(res_indices) == 0: continue

        Try = min([20, len(res_indices)])
        minDist, targetX, targetY, dist, ct = 100, 0, 0, 0, 0

        while ct < Try and RIDs is not []:
            ct += 1
            r = len(res_indices)
            Try = min([20, r])
            j = randint(0, r-1)

            if fxn == 'chemotaxis':
                Rtype = Rtypes[j]
                x2 = RX[j]
                y2 = RY[j]
            elif fxn == 'pred-prey':
                Rtype = SpeciesIDs[j]
                x2 = IndX[j]
                y2 = IndY[j]

            if restype == Rtype: # is capable of consuming the resource type

                dist = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

                if dist < minDist:
                    minDist = dist
                    if fxn == 'chemotaxis':
                        targetX = RX[j]
                        targetY = RY[j]
                    elif fxn == 'pred-prey':
                        targetX = IndX[j]
                        targetY = IndY[j]

        # A cost for active dispersal
        Q = CellQuotas[i]
        maint = MaintDict[spID]
        mfd = MainFactorDict[spID]

        Q -= maint*disp

        if Q < maint/mfd:
            IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
            ResLists = Rvals, Rtypes, RX, RY, RIDs
            IndLists, ResLists, RID = dead(IndLists, ResLists, RID, i, Q, TC, RC)
            numDead += 1
            continue

        if x1 > targetX:
            x1 -= dist*disp
        elif x1 < targetX:
            x1 += dist*disp
        if y1 > targetY:
            y1 -= dist*disp
        elif y1 < targetY:
            y1 += dist*disp

        IndX[i] = x1
        IndY[i] = y1
        CellQuotas[i] = Q

    IndLists = SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList
    ResLists = Rvals, Rtypes, RX, RY, RIDs
    return [IndLists, IndDicts, IndID, ResLists, ResDict, RID, numDead]
