# -*- coding: utf-8 -*-
from __future__ import division
from random import randint, choice
import numpy as np
import math
import time
import sys
import os

mydir = os.path.expanduser("~/GitHub/Micro-Encounter")

def checkVal(val, line):
    if val < 0:
        print 'line',line,': error: val < 0:', val
        sys.exit()
    return


def get_closest(RIDs, RX, RY, RZ, Rtypes, restype, coords):

    closest = 'none'
    res_indices = [xi for xi, x in enumerate(Rtypes) if x == restype or x == 'dead']

    if len(res_indices) == 0:
        return closest

    x1, y1, z1 = coords
    Try = min([40, len(res_indices)])
    minDist, ct = 10**10, 0

    while ct < Try:
        ct += 1
        j = randint(0, len(res_indices)-1)
        x = RX[j]
        y = RY[j]
        z = RZ[j]

        dist = math.sqrt((x1 - x)**2 + (y1 - y)**2 + (z1 - z)**2)

        if dist < minDist:
            minDist = dist
            closest = RIDs[j]

        return closest



def dead(IndDict, ResLists, i, Q, TC, RC):
    Rvals, Rtypes, RX, RY, RZ, RIDs = ResLists

    Q = IndDict[i]['quota']
    x = IndDict[i]['x']
    y = IndDict[i]['y']
    z = IndDict[i]['z']

    if '-scavenging-' in TC:
        Rvals.append(Q)
        Rtypes.append('dead')
        RX.append(x)
        RY.append(y)
        RZ.append(z)
        RIDs.append(time.clock())

    del IndDict[i]

    ResLists = Rvals, Rtypes, RX, RY, RZ, RIDs
    return [IndDict, ResLists]



def ResIn(ResLists, ResDict, params, ct, ComplexityLevels):

    Rvals, Rtypes, RX, RY, RZ, RIDs = ResLists
    width, height, length, seed, m, r, gmax, mmax, dmax, pmax, mfact, std = params
    SC, TC, RC = ComplexityLevels

    Ymean = np.random.uniform(0, height)
    Xmean = np.random.uniform(0, width)
    Zmean = np.random.uniform(0, length)
    sizemax = 10000

    res_in = 0
    x = np.random.binomial(1, r)
    if x == 1:
        res_in += 1
        RIDs.append(time.clock())
        Rvals.append(randint(sizemax, sizemax))
        RY.append(np.random.uniform(0, height)) # default is random
        RX.append(np.random.uniform(0, width)) # default is random
        RZ.append(np.random.uniform(0, length)) # default is random
        pmin, pmax = 1.0, 1.0 # default is 'simple'
        ri = 'a' # default is monoculture

        if '-polyculture-' in RC:
            ri = choice(['a', 'b', 'c']) #ri = choice(string.lowercase) # default is polyculture

        if '-lockandkey-' in RC:
            pmin, pmax = 0.1, 0.5

        RY[-1] = np.random.normal(Ymean, std)
        RX[-1] = np.random.normal(Xmean, std)
        RZ[-1] = np.random.normal(Zmean, std)

        if ri not in ResDict:
            ResDict[ri] = np.random.uniform(pmin, pmax)
        Rtypes.append(ri)

    ResLists = Rvals, Rtypes, RX, RY, RZ, RIDs
    return [ResLists, ResDict, res_in]



def immigration(IndDict, SpDicts, params, ct, ComplexityLevels):

    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = SpDicts
    width, height, length, seed, m, r, gmax, mmax, dmax, pmax, mfact, std = params
    SC, TC, RC = ComplexityLevels

    Ymean = np.random.uniform(0, height)
    Xmean = np.random.uniform(0, width)
    Zmean = np.random.uniform(0, length)

    if ct > 1:
        seed = np.random.binomial(1, m)

    for num in range(seed):
        IndID = time.clock()
        spID = np.random.randint(1, 100)
        IndDict[IndID] = {'species' : int(spID)}
        IndDict[IndID]['quota'] = np.random.uniform(100, 1000)
        IndDict[IndID]['state'] = choice(['active'])
        IndDict[IndID]['closest'] = 'none'
        IndDict[IndID]['xdir'] = choice([1, -1])
        IndDict[IndID]['ydir'] = choice([1, -1])
        IndDict[IndID]['zdir'] = choice([1, -1])

        vals = np.random.normal([Ymean, Xmean, Zmean], std)

        IndDict[IndID]['y'] = vals[0]
        IndDict[IndID]['x'] = vals[1]
        IndDict[IndID]['z'] = vals[2]

        if spID not in GrowthDict:
            TrophicDict[spID] = {'resource': str(), 'crossfeed_to': str()}

            if '-crossfeeding-' in TC:
                rlist = ['a', 'b', 'c']
                r = choice(rlist)
                rlist.pop(rlist.index(r))

                TrophicDict[spID]['resource'] = r
                #TrophicDict[prop]['resource'] = choice(string.lowercase)

                TrophicDict[spID]['crossfeeds_to'] = choice(rlist)
                #TrophicDict[prop]['crossfeeds_to'] = choice(string.lowercase)

            elif '-none-' in TC or '-scavenging-' in TC:
                if '-monoculture-' in RC:
                    TrophicDict[spID]['resource'] = 'a'
                elif '-polyculture-' in RC:
                    TrophicDict[spID]['resource'] = choice(['a', 'b', 'c'])

            # species growth rate
            GrowthDict[spID] = np.random.uniform(gmax/10, gmax)

            # species maintenance
            MaintDict[spID] = np.random.uniform(mmax/10, mmax)

            # species maintenance reduction factor
            MainFactorDict[spID] = np.random.uniform(mfact/10, mfact)

            # species resuscitation factor
            RPFDict[spID] = np.random.uniform(pmax/10, pmax)

            # species active dispersal rate
            DispDict[spID] = np.random.uniform(dmax/10, dmax)

    SpDicts = GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict
    return [IndDict, SpDicts]




def maintenance(IndDict, SpDicts, ResLists, ResDict, ComplexityLevels, numDead):

    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = SpDicts
    Rvals, Rtypes, RX, RY, RZ, RIDs = ResLists
    SC, TC, RC = ComplexityLevels

    IndIDs = list(IndDict.keys())
    n = len(IndIDs)
    if n == 0:
        return [IndDict, ResLists, ResDict, numDead]

    while len(IndIDs) > 0:
        i = choice(IndIDs)
        Q = IndDict[i]['quota']
        state = IndDict[i]['state']
        spID = IndDict[i]['species']
        mfd = MainFactorDict[spID]
        maint = float(MaintDict[spID])

        if state == 'dormant':
            maint = maint/mfd

        if Q > maint:
            Q -= maint
            IndDict[i]['quota'] = Q

        else:
            ResLists = Rvals, Rtypes, RX, RY, RZ, RIDs
            IndDict, ResLists = dead(IndDict, ResLists, i, Q, TC, RC)
            numDead += 1

        IndIDs.pop(IndIDs.index(i))

    ResLists = Rvals, Rtypes, RX, RY, RZ, RIDs
    SpDicts = GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict
    return [IndDict, ResLists, ResDict, numDead]



def transition(IndDict, SpDicts, ResLists, ResDict, ComplexityLevels, numDead):

    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = SpDicts
    Rvals, Rtypes, RX, RY, RZ, RIDs = ResLists
    SC, TC, RC = ComplexityLevels

    IndIDs = list(IndDict.keys())
    n = len(IndIDs)
    if n == 0:
        return [IndDict, ResLists, ResDict, numDead]

    while len(IndIDs) > 0:
        i = choice(IndIDs)
        IndIDs.pop(IndIDs.index(i))

        Q = IndDict[i]['quota']
        state = IndDict[i]['state']
        spID = IndDict[i]['species']
        maint = float(MaintDict[spID])
        rpd = RPFDict[spID]

        x = np.random.binomial(1, rpd)
        if state == 'dormant' and x == 1:
            IndDict[i]['state'] = 'active'

        elif state == 'active' and Q <= maint*4:
            IndDict[i]['state'] = 'dormant'

        else: continue

        if Q < maint:
            ResLists = Rvals, Rtypes, RX, RY, RZ, RIDs
            IndDict, ResLists = dead(IndDict, ResLists, i, Q, TC, RC)
            numDead += 1

    ResLists = Rvals, Rtypes, RX, RY, RZ, RIDs
    return [IndDict, ResLists, ResDict, numDead]




def consume(IndDict, SpDicts, ResLists, ResDict, params, ComplexityLevels, numDead):

    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = SpDicts
    width, height, length, seed, m, r, gmax, mmax, dmax, pmax, mfact, std = params
    Rvals, Rtypes, RX, RY, RZ, RIDs = ResLists
    SC, TC, RC = ComplexityLevels

    maxQ, encounters = 1000, 0
    IndIDs = IndDict.keys()
    n = len(IndIDs)
    r = len(RIDs)
    if n >= 0 and r == 0:
        return [ResLists, ResDict, IndDict, encounters, numDead]

    while len(IndIDs) > 0:
        i = choice(IndIDs)
        IndIDs.pop(IndIDs.index(i))
        Q = IndDict[i]['quota']
        spID = IndDict[i]['species']
        state = IndDict[i]['state']
        x1 = IndDict[i]['x']
        y1 = IndDict[i]['y']
        z1 = IndDict[i]['z']

        closest = IndDict[i]['closest']

        if closest not in RIDs or 'randwalk' in SC:
            restype = TrophicDict[spID]['resource']

            coords = [x1, y1, z1]
            IndDict[i]['closest'] = get_closest(RIDs, RX, RY, RZ, Rtypes, restype, coords)
            closest = IndDict[i]['closest']
            if closest == 'none':
                continue

        j = RIDs.index(closest)
        Rval = Rvals[j]
        Rtype = Rtypes[j]
        x2 = RX[j]
        y2 = RY[j]
        z2 = RZ[j]

        dist = math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)

        i_radius = ((0.75*Q)/math.pi)**(1.0/3)
        r_radius = 100 * ((0.75*Rval)/math.pi)**(1.0/3)

        #print 'dists', dist, r_radius, i_radius, 'Rtype:', Rtype, SC
        if dist <= i_radius + r_radius:
            if state == 'dormant':
                IndDict[i]['state'] = 'active'

            #print 'dists', dist, r_radius, i_radius, 'Rtype:', Rtype, SC
            if np.random.binomial(1, ResDict[Rtype]) == 1:
                encounters += 1
                mu = GrowthDict[spID] * Q

                if mu >= Rval:
                    Q += Rval
                    Rval = 0.0

                elif mu < Rval:
                    Q += mu
                    Rval -= mu

                if Q > maxQ:
                    Rval += Q - maxQ
                    Q = float(maxQ)

                IndDict[i]['quota'] = Q
                checkVal(Q, 348)

                if Rval == 0:
                    Rvals.pop(j)
                    RIDs.pop(j)
                    Rtypes.pop(j)
                    RX.pop(j)
                    RY.pop(j)
                    RZ.pop(j)

                    if len(RIDs) == 0:
                        ResLists = Rvals, Rtypes, RX, RY, RZ, RIDs
                        return [ResLists, ResDict, IndDict, encounters, numDead]

                else:
                    r = np.random.uniform(0, Rval)
                    rlist = [r, Rval - r]

                    Rvals[j] = max(rlist)
                    r2 = min(rlist)

                    if '-crossfeeding-' in TC:
                        Rtype = TrophicDict[spID]['crossfeeds_to']
                        pmin, pmax = 1.0, 1.0 # 'simple' is default

                        if '-lockandkey-' in RC:
                            pmin, pmax = 0.1, 0.5

                        if Rtype not in ResDict:
                            ResDict[Rtype] = np.random.uniform(pmin, pmax)

                    RIDs.append(time.clock())
                    Rtypes.append(Rtype)
                    Rvals.append(r2)
                    RX.append(x2)
                    RY.append(y2)
                    RZ.append(z2)


    ResLists = Rvals, Rtypes, RX, RY, RZ, RIDs
    return [ResLists, ResDict, IndDict, encounters, numDead]



def reproduce(IndDict, SpDicts, ResLists, ResDict, params, ComplexityLevels, numDead):

    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = SpDicts
    width, height, length, seed, m, r, gmax, mmax, dmax, pmax, mfact, std = params
    Rvals, Rtypes, RX, RY, RZ, RIDs = ResLists
    SC, TC, RC = ComplexityLevels

    maxQ, PRODI = 1000, 0
    IndIDs = list(IndDict.keys())
    n = len(IndIDs)
    if n == 0: return [PRODI, IndDict, ResLists, ResDict, numDead]

    while len(IndIDs) > 0:
        i = choice(IndIDs)
        IndIDs.pop(IndIDs.index(i))

        age = time.clock() - i
        Q = IndDict[i]['quota']
        spID = IndDict[i]['species']
        state = IndDict[i]['state']
        i_radius = ((0.75*Q)/math.pi)**(1.0/3)

        if state == 'dormant' or age < 0.12: continue

        elif i_radius > 5.0 and np.random.binomial(1, Q/maxQ) == 1: # individual is large enough to reproduce
            PRODI += 1
            i2 = time.clock()

            IndDict[i]['quota'] = Q/2
            IndDict[i2] = {'quota': Q/2}
            IndDict[i2]['species'] = spID
            IndDict[i2]['state'] = 'active'
            IndDict[i2]['xdir'] = choice([1, -1])
            IndDict[i2]['ydir'] = choice([1, -1])
            IndDict[i2]['zdir'] = choice([1, -1])
            IndDict[i2]['x'] = IndDict[i]['x']
            IndDict[i2]['y'] = IndDict[i]['y']
            IndDict[i2]['z'] = IndDict[i]['z']
            IndDict[i2]['closest'] = IndDict[i]['closest']

    ResLists = Rvals, Rtypes, RX, RY, RZ, RIDs
    return [PRODI, IndDict, ResLists, ResDict, numDead]




def res_dispersal(ResLists, params, ct, ComplexityLevels):

    Rvals, Rtypes, RX, RY, RZ, RIDs = ResLists
    width, height, length, seed, m, r, gmax, mmax, dmax, pmax, mfact, std = params
    SC, TC, RC = ComplexityLevels

    n = len(RIDs)
    if n == 0: return [Rvals, Rtypes, RX, RY, RZ, RIDs]

    for j in range(n):
        i = randint(0, len(RIDs)-1)

        x = RX[i]
        y = RY[i]
        z = RZ[i]
        rval = Rvals[i]

        if rval <= 0:
            Rvals.pop(i)
            RIDs.pop(i)
            Rtypes.pop(i)
            RX.pop(i)
            RY.pop(i)
            RZ.pop(i)
            continue

        n = len(RIDs)
        if n == 0: return [Rvals, Rtypes, RX, RY, RZ, RIDs]

        if '-wellmixed-' in SC:
            y = float(np.random.uniform(0, height))
            x = float(np.random.uniform(0, width))
            z = float(np.random.uniform(0, length))

        elif '-brownian-' in SC:
            y += choice([-1,1]) * float(np.random.uniform(0, 1))
            x += choice([-1,1]) * float(np.random.uniform(0, 1))
            z += choice([-1,1]) * float(np.random.uniform(0, 1))

            if y < 0: y = 0
            elif y > height: y = height
            if x < 0: x = 0
            elif x > width: x = width
            if z < 0: z = 0
            elif z > length: z = 0

        RX[i] = x
        RZ[i] = z
        RY[i] = y

    ResLists = Rvals, Rtypes, RX, RY, RZ, RIDs
    return ResLists



def dispersal(IndDict, SpDicts, ResLists, ResDict, params, ComplexityLevels, numDead):

    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = SpDicts
    width, height, length, seed, m, r, gmax, mmax, dmax, pmax, mfact, std = params
    Rvals, Rtypes, RX, RY, RZ, RIDs = ResLists
    SC, TC, RC = ComplexityLevels

    n = len(IndDict.keys())
    if n == 0: return [IndDict, ResLists, ResDict, numDead]

    if '-chemotaxis-' in SC:
        IndDict, ResLists, ResDict, numDead = chemotaxis(IndDict, SpDicts, ResLists, ResDict, params, ComplexityLevels, numDead)

    elif '-wellmixed-' in SC or '-randwalk-' in SC:
        IndIDs = IndDict.keys()

        while len(IndIDs) > 0:

            i = choice(IndIDs)
            IndIDs.pop(IndIDs.index(i))

            Q = IndDict[i]['quota']
            spID = IndDict[i]['species']
            x1 = IndDict[i]['x']
            y1 = IndDict[i]['y']
            z1 = IndDict[i]['z']
            xdir = IndDict[i]['xdir']
            ydir = IndDict[i]['ydir']
            zdir = IndDict[i]['zdir']

            maint = MaintDict[spID]
            disp = DispDict[spID]
            dist = 0.0

            if '-wellmixed-' in SC:
                x = np.random.uniform(0, width)
                y = np.random.uniform(0, height)
                z = np.random.uniform(0, length)
                continue

            elif '-randwalk-' in SC:
                x = x1 + xdir * np.random.uniform(0, disp*width)
                y = y1 + ydir * np.random.uniform(0, disp*height)
                z = z1 + zdir * np.random.uniform(0, disp*length)

                dist = math.sqrt((x1 - x)**2 + (y1 - y)**2 + (z1 - z)**2)

            if x > width:
                x = width
                xdir = -1*xdir
            elif x < 0:
                x = 0
                xdir = -1*xdir

            if y > height:
                y = height
                ydir = -1*ydir
            elif y < 0:
                y = 0
                ydir = -1*ydir

            if z > length:
                z = length
                zdir = -1*zdir
            elif z < 0:
                z = 0
                zdir = -1*zdir

            if '-randwalk-' in SC:
                # A cost for active dispersal
                if Q - maint*(dist/width) <= 0.0:
                    ResLists = Rvals, Rtypes, RX, RY, RZ, RIDs
                    IndDict, ResLists = dead(IndDict, ResLists, i, Q, TC, RC)
                    numDead += 1
                else:
                    Q -= maint*(dist/width)
                    checkVal(Q, 470)
                    IndDict[i]['quota'] = Q
                    IndDict[i]['x'] = x
                    IndDict[i]['xdir'] = xdir
                    IndDict[i]['y'] = y
                    IndDict[i]['ydir'] = ydir
                    IndDict[i]['z'] = z
                    IndDict[i]['zdir'] = zdir


    ResLists = Rvals, Rtypes, RX, RY, RZ, RIDs
    return [IndDict, ResLists, ResDict, numDead]




def chemotaxis(IndDict, SpDicts, ResLists, ResDict, params, ComplexityLevels, numDead):

    GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = SpDicts
    width, height, length, seed, m, r, gmax, mmax, dmax, pmax, mfact, std = params
    Rvals, Rtypes, RX, RY, RZ, RIDs = ResLists
    SC, TC, RC = ComplexityLevels

    IndIDs = list(IndDict.keys())
    n = len(IndIDs)
    if n == 0: return [IndDict, ResLists, ResDict, numDead]

    while len(IndIDs) > 0:

        i = choice(IndIDs)
        IndIDs.pop(IndIDs.index(i))
        state = IndDict[i]['state']
        if state == 'dormant': continue

        Q = IndDict[i]['quota']
        x1 = IndDict[i]['x']
        y1 = IndDict[i]['y']
        z1 = IndDict[i]['z']
        xdir = IndDict[i]['xdir']
        ydir = IndDict[i]['ydir']
        zdir = IndDict[i]['zdir']

        spID = IndDict[i]['species']
        maint = MaintDict[spID]
        disp = DispDict[spID]

        restype = TrophicDict[spID]['resource']
        closest = IndDict[i]['closest']

        if closest not in RIDs:
            coords = [x1, y1, z1]
            IndDict[i]['closest'] = get_closest(RIDs, RX, RY, RZ, Rtypes, restype, coords)

        if closest not in RIDs:
            x = x1 + xdir * np.random.uniform(0, disp*width)
            y = y1 + ydir * np.random.uniform(0, disp*height)
            z = z1 + zdir * np.random.uniform(0, disp*length)

            dist = math.sqrt((x1 - x)**2 + (y1 - y)**2 + (z1 - z)**2)

            if x > width:
                x = width
                xdir = -1*xdir
            elif x < 0:
                x = 0
                xdir = -1*xdir

            if y > height:
                y = height
                ydir = -1*ydir
            elif y < 0:
                y = 0
                ydir = -1*ydir

            if z > length:
                z = length
                zdir = -1*zdir
            elif z < 0:
                z = 0
                zdir = -1*zdir

            # A cost for active dispersal
            if Q - maint*(dist/width) <= 0.0:
                ResLists = Rvals, Rtypes, RX, RY, RZ, RIDs
                IndDict, ResLists = dead(IndDict, ResLists, i, Q, TC, RC)
                numDead += 1
            else:
                Q -= maint*(dist/width)
                checkVal(Q, 470)
                IndDict[i]['quota'] = Q
                IndDict[i]['x'] = x
                IndDict[i]['xdir'] = xdir
                IndDict[i]['y'] = y
                IndDict[i]['ydir'] = ydir
                IndDict[i]['z'] = z
                IndDict[i]['zdir'] = zdir

        elif closest in RIDs:
            ri = RIDs.index(closest)
            x2 = RX[ri]
            y2 = RY[ri]
            z2 = RZ[ri]

            x = np.abs(x1 - x2)
            if x1 > x2:
                x1 -= np.random.uniform(0, disp*x)
            elif x1 < x2:
                x1 += np.random.uniform(0, disp*x)

            y = np.abs(y1 - y2)
            if y1 > y2:
                y1 -= np.random.uniform(0, disp*y)
            elif y1 < y2:
                y1 += np.random.uniform(0, disp*y)

            z = np.abs(z1 - z2)
            if z1 > z2:
                z1 -= np.random.uniform(0, disp*z)
            elif z1 < z2:
                z1 += np.random.uniform(0, disp*z)

            dist = math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)

            if Q - maint*(dist/width)*1.0 <= 0.0:
                ResLists = Rvals, Rtypes, RX, RY, RZ, RIDs
                IndDict, ResLists = dead(IndDict, ResLists, i, Q, TC, RC)
                numDead += 1

            else:
                Q -= maint*(dist/width)*1.0
                checkVal(Q, 657)
                IndDict[i]['quota'] = Q
                IndDict[i]['x'] = x1
                IndDict[i]['y'] = y1
                IndDict[i]['z'] = z1

    ResLists = Rvals, Rtypes, RX, RY, RZ, RIDs
    return [IndDict, ResLists, ResDict, numDead]
