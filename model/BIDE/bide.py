# -*- coding: utf-8 -*-
from __future__ import division
from random import randint, choice
import numpy as np
import sys
import math
import os

mydir = os.path.expanduser("~/")
sys.path.append(mydir + "GitHub/Micro-Encounter/model/EnzymeField")
import EnzymeField as field


limit = 0.1


def coord(d):
    return float(np.random.uniform(0.1*d, 0.9*d))


def GetIndParam(means):
    vals = []

    if isinstance(means, float) or isinstance(means, int):
        std = means/1000.0
        vals = np.random.normal(means, std)
        if vals < 0.00001:
            vals = 0.00001

    else:
        for val in means:
            std = val/1000.0
            i = np.random.normal(val, std)
            if i < 0.00001:
                i = 0.00001
            vals.append(i)

    return vals



def GetRAD(vector):
    RAD = []
    unique = list(set(vector))

    for val in unique:
        RAD.append(vector.count(val)) # the abundance of each Species

    return RAD, unique # the rad and the specieslist


def ResIn(std, ct, RList, Vals, Xs, Ys, ID, IDs, numr, w, h, u0, TrophicComplexityLevel, SpatialComplexityLevel, ResourceComplexityLevel, BiologicalComplexityLevel):

    Ymean = float(np.random.uniform(0.001*h, 0.999*h))
    Xmean = float(np.random.uniform(0.001*w, 0.999*w))

    listlen = [len(RList), len(Vals), len(Xs), len(Ys), len(IDs)]
    if min(listlen) != max(listlen):
        print ct, 'In ResIn (top)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()

    for r in range(numr):
        x = np.random.binomial(1, u0)

        if x == 1:
            IDs.append(ID)
            ID += 1

            if ResourceComplexityLevel == 1:
                #http://stackoverflow.com/questions/1957273/how-do-i-generate-a-random-string-of-length-x-a-z-only-in-python
                # all resource particles are of the same type
                r = choice(['a', 'aa', 'aaaa'])
                if r == 'a':
                    n = randint(1, 40)
                elif r == 'aa':
                    n = randint(1, 20)
                elif r == 'aaaa':
                    n = randint(1, 10)
                r = ''.join(r for _ in xrange(n))

            elif ResourceComplexityLevel == 2:
                # Three types of resource particles.
                r = choice(['a', 'bb', 'cccc'])
                if r == 'a':
                    n = randint(1, 40)
                elif r == 'bb':
                    n = randint(1, 20)
                elif r == 'cccc':
                    n = randint(1, 10)
                r = ''.join(r for _ in xrange(n))

            elif ResourceComplexityLevel == 3:
                # Three types of resource particles that vary in the effort/energy needed to break them down
                r = choice(['a-', 'b-b-', 'c-c-c-c-'])
                if r == 'a-':
                    n = randint(1, 40)
                elif r == 'b-b-':
                    n = randint(1, 20)
                elif r == 'c-c-c-c-':
                    n = randint(1, 10)

                r = ''.join(r for _ in xrange(n))

            RList.append(r)
            Rtype = r[0]
            Vals.append(r.count(Rtype))

            if SpatialComplexityLevel == 1:

                Ys.append(float(np.random.uniform(0.001*h, 0.999*h)))
                Xs.append(float(np.random.uniform(0.001*w, 0.999*w)))

            elif SpatialComplexityLevel <= 3:

                vals = np.random.normal([Ymean, Xmean], std)

                if vals[0] < 0.001*h:
                    vals[0] = 0.001*h
                elif vals[0] > 0.999*h:
                    vals[0] = 0.999*h

                if vals[1] < 0.001*w:
                    vals[1] = 0.001*w
                elif vals[1] > 0.999*w:
                    vals[1] = 0.999*w

                Ys.append(vals[0])
                Xs.append(vals[1])

    listlen = [len(RList), len(Vals), len(Xs), len(Ys), len(IDs)]
    if min(listlen) != max(listlen):
        print 'In ResIn (bottom)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()

    return [RList, Vals, Xs, Ys, IDs, ID]



def immigration(std, mfmax, p_max, d_max, g_max, m_max, seed, ip, Sp, Xs, Ys, w, h, MD,
    MFD, RPD, EnvD, envGs, GD, DispD, IDs, ID, Qs, RD, u0, alpha, GList, MList,
    DList, ADList, EVList, TLList, ct, TrophicComplexityLevel, SpatialComplexityLevel,
    ResourceComplexityLevel, BiologicalComplexityLevel):

    Ymean = float(np.random.uniform(0.001*h, 0.999*h))
    Xmean = float(np.random.uniform(0.001*w, 0.999*w))

    listlen = [len(Sp), len(Qs), len(IDs), len(Xs), len(Ys), len(GList), len(MList), len(DList), len(ADList), len(EVList), len(TLList)]
    if min(listlen) != max(listlen):
        print ct, 'In immigration (top)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()


    if ct > 1: seed = 1
    if ct == 1: SpatialComplexityLevel = 1

    for m in range(seed):
        x = 0

        if seed > 1:
            x = 1

        else:
            x = np.random.binomial(1, ip)

        if x == 1:
            prop = str(float(np.random.logseries(alpha, 1)))

            Sp.append(prop)

            if SpatialComplexityLevel <= 3:
                Ys.append(float(np.random.uniform(0.001*h, 0.999*h)))
                Xs.append(float(np.random.uniform(0.001*w, 0.999*w)))

            elif SpatialComplexityLevel <= 3:

                vals = np.random.normal([Ymean, Xmean], std)

                if vals[0] < 0.001*h:
                    vals[0] = 0.001*h
                elif vals[0] > 0.999*h:
                    vals[0] = 0.999*h

                if vals[1] < 0.001*w:
                    vals[1] = 0.001*w
                elif vals[1] > 0.999*w:
                    vals[1] = 0.999*w

                Ys.append(vals[0])
                Xs.append(vals[1])

            IDs.append(ID)
            ID += 1

            tl = choice(['a', 'b', 'c'])
            TLList.append(tl)

            Q = float(np.random.uniform(0.01, 0.1))
            Qs.append(Q)

            if prop not in GD:

                # species growth rate
                #g = np.random.uniform(g_max/1, g_max)
                GD[prop] = g_max

                # species maintenance
                MD[prop] = m_max #np.random.uniform(m_max/1, m_max)

                # species maintenance reduction factor
                MFD[prop] = mfmax #np.random.uniform(1, mfmax)

                # species resuscitation factor
                RPD[prop] = p_max #np.random.uniform(0.1, 0.1)

                # species active dispersal rate
                DispD[prop] = d_max #np.random.uniform(d_max/1, d_max)

                # species environmental gradient optima
                glist = []
                for e in envGs:
                    x = np.random.uniform(0.0, w)
                    y = np.random.uniform(0.0, h)
                    glist.append([x,y])
                EnvD[prop] = glist

                # Resource use efficiency
                RD[prop] = [1.0 - g_max]

            ADList.append('a')

            means = GD[prop]
            i = GetIndParam(means)
            GList.append(i)

            means = MD[prop]
            i = GetIndParam(means)
            MList.append(i)

            means = RD[prop]
            n = GetIndParam(means)
            EVList.append(n) # resource efficiency value list

            means = DispD[prop]
            i = GetIndParam(means)
            DList.append(i)


    listlen = [len(Sp), len(Qs), len(IDs), len(Xs), len(Ys), len(GList), len(MList), len(DList), len(ADList), len(EVList), len(TLList)]
    if min(listlen) != max(listlen):
        print 'In immigration (bottom)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()


    return [Sp, Xs, Ys, MD, MFD, RPD, EnvD, GD, DispD, IDs, ID, Qs,
            RD, GList, MList, DList, ADList, EVList, TLList]






def maintenance(SpeciesIDs, Xs, Ys, MD, MFD, RPD, EnvD, IDs, Qs, GList, MaintList,
    DList, ADList, EVList, TLList, RList, RVals, RXs, RYs, RIDs, RID,
    TrophicComplexityLevel, SpatialComplexityLevel,ResourceComplexityLevel, BiologicalComplexityLevel):

    listlen = [len(SpeciesIDs), len(Qs), len(IDs), len(Xs), len(Ys), len(GList), len(MaintList), len(DList), len(ADList), len(EVList), len(TLList)]
    if min(listlen) != max(listlen):
        print 'In maintenance (top)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()

    if SpeciesIDs == []:
        return [SpeciesIDs, Xs, Ys, IDs, Qs, GList, MaintList, DList, ADList, EVList, TLList, RList, RVals, RXs, RYs, RIDs, RID]

    n = len(IDs)
    for j in range(n):

        i = randint(0, len(IDs)-1)

        val = Qs[i]
        val -= MaintList[i] # maintanence influenced by species id

        if val <= MaintList[i]*0.001:   # starved

            if TrophicComplexityLevel >= 3:
                RList.append('d')
                RVals.append(1)
                RXs.append(Xs[i])
                RYs.append(Ys[i])
                RIDs.append(RID)
                RID += 1

            Qs.pop(i)
            SpeciesIDs.pop(i)
            IDs.pop(i)
            Xs.pop(i)
            Ys.pop(i)
            GList.pop(i)
            MaintList.pop(i)
            TLList.pop(i)
            DList.pop(i)
            ADList.pop(i)
            EVList.pop(i)

        else: Qs[i] = val

    listlen = [len(SpeciesIDs), len(Qs), len(IDs), len(Xs), len(Ys), len(GList), len(MaintList), len(DList), len(ADList), len(EVList), len(TLList)]
    if min(listlen) != max(listlen):
        print 'In maintenance (bottom)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()

    return [SpeciesIDs, Xs, Ys, IDs, Qs, GList, MaintList, DList, ADList, EVList, TLList, RList, RVals, RXs, RYs, RIDs, RID]





def transition(SpeciesIDs, Xs, Ys, GList, DList, ADList, EVList, IDs, Qs, MaintList, TLList, RList, RVals, RXs, RYs, RIDs, RID, MFD, RPD, TrophicComplexityLevel, SpatialComplexityLevel, ResourceComplexityLevel, BiologicalComplexityLevel):

    if SpeciesIDs == []:
        return [SpeciesIDs, Xs, Ys, GList, DList, ADList, EVList, IDs, Qs, GList, MaintList, TLList, RList, RVals, RXs, RYs, RIDs, RID]

    n = len(IDs)
    for j in range(n):

        i = randint(0, len(IDs)-1)
        spid = SpeciesIDs[i]
        state = ADList[i]

        mfd = MFD[spid]  # multiplicative maintenance factor (is greater than 1)

        if state == 'd':

            x = np.random.binomial(1, RPD[spid]) # make this probability a randomly chosen variable
            if x == 1:
                ADList[i] = 'a'
                MaintList[i] = MaintList[i]*mfd # metabolic maintenance increases

        if state == 'a':

            Q = Qs[i]
            if Q <= MaintList[i]*mfd:  # go dormant
                MaintList[i] = MaintList[i]/mfd # metabolic maintenance decreases
                ADList[i] = 'd'

            if Q < 0.0:

                Qs.pop(i)
                SpeciesIDs.pop(i)
                IDs.pop(i)
            	Xs.pop(i)
                Ys.pop(i)
                GList.pop(i)
                MaintList.pop(i)
                TLList.pop(i)
                DList.pop(i)
            	ADList.pop(i)
                EVList.pop(i)
                continue


    return [SpeciesIDs, Xs, Ys, GList, DList, ADList, EVList, IDs, Qs, GList, MaintList, TLList, RList, RVals, RXs, RYs, RIDs, RID]





def consume(field, RList, RVals, RIDs, RID, RXs, RYs, SpeciesIDs, Qs, IndIDs, IndID,
    IXs, IYs, w, h, GD, RD, DispD, GrowthList, MaintList, DispList, ADList, TLList, EVList,
    TrophicComplexityLevel, SpatialComplexityLevel, ResourceComplexityLevel, BiologicalComplexityLevel, enz_action = True):

    encounters = 0

    listlen = [len(SpeciesIDs), len(Qs), len(IndIDs), len(IXs), len(IYs), len(GrowthList), len(MaintList), len(DispList), len(ADList), len(EVList), len(TLList)]

    if min(listlen) != max(listlen):
        print 'In consume (top)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()

    if SpeciesIDs == [] or RIDs == []:
        return [RList, RVals, RIDs, RID, RXs, RYs, SpeciesIDs, Qs, encounters]

    n = len(IndIDs)

    for ii in range(n):
        i = randint(0, n-1)

        # Trophic level
        tl = TLList[i]
        Q = Qs[i]

        x1 = IXs[i]
        y1 = IYs[i]

        Try = min([50, len(RIDs)])
        ct = 0

        while ct < Try and RIDs != []:
            ct += 1

            r = len(RIDs)
            Try = min([50, r])

            j = randint(0, r-1)
            # The food
            R = RList[j]

            if R.startswith('-'):
                R = R[1:]
            if R.endswith('-'):
                R = R[:-1]

            Rtype = R[0]

            if Rtype == '-':
                print 'Rtype is a hyphen'
                sys.exit()

            if ResourceComplexityLevel == 1 or Rtype == 'd':
                tl = str(Rtype)

            if tl == Rtype: # individual is capable of consuming the resource type

                x2 = RXs[j]
                y2 = RYs[j]
                dist = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

                ind_radius = np.mean(Qs[i])
                Rval = RVals[j]
                res_radius = Rval/20.0

                if dist <= ind_radius + res_radius:
                    ADList[i] = 'a'

                    #encounters += 1
                    eat = 'no'
                    ct = Try

                    if '-' in R:
                        # Useful stackoverflow answers:
                        # http://stackoverflow.com/questions/32212998/split-string-by-n-when-n-is-random
                        # http://stackoverflow.com/questions/19954593/python-checking-a-strings-first-and-last-character
                        # http://stackoverflow.com/questions/5188792/how-to-check-a-string-for-specific-characters

                        pos = randint(1, len(R)-1)

                        if R[pos] == '-':
                            pieces = [R[:pos], R[pos:]]

                            for index, piece in enumerate(pieces):
                                if piece.startswith('-'):
                                    piece = piece[1:]

                                if piece.endswith('-'):
                                    piece = piece[:-1]

                            if '-' in pieces[0] and '-' in pieces[1]:

                                RList[j] = pieces[0]
                                RVals[j] = pieces[0].count(Rtype)

                                RList.append(pieces[1])
                                RVals.append(pieces[1].count(Rtype))
                                RIDs.append(RID)
                                RID += 1
                                RXs.append(x2)
                                RYs.append(y2)

                            else:
                                for index, piece in enumerate(pieces):

                                    if '-' not in piece:
                                        eat = 'yes'

                                        RList[j] = piece
                                        RVals[j] = piece.count(Rtype)

                                        R = pieces[index]
                                        Rval = pieces[index].count(Rtype)

                                        RList.append(pieces[index-1])
                                        RVals.append(pieces[index-1].count(Rtype))
                                        RIDs.append(RID)
                                        RID += 1
                                        RXs.append(x2)
                                        RYs.append(y2)
                                        break


                    elif '-' not in R and ResourceComplexityLevel < 3:
                        eat = 'yes'

                        if len(R) > 1:

                            pieces = [R[:1], R[1:]] # individual consumes one letter at a time

                            RList[j] = pieces[0]
                            RVals[j] = len(pieces[0])

                            RList.append(pieces[1])
                            RVals.append(len(pieces[1]))
                            RIDs.append(RID)
                            RID += 1
                            RXs.append(x2)
                            RYs.append(y2)

                    if eat == 'yes':

                        encounters += 1
                        # The individual's cell quota
                        Q = Qs[i]

                        # The species
                        sp = SpeciesIDs[i]
                        mu = GD[sp]

                        # Increase cell quota & decrease resource particle size
                        Q += (mu * Q)
                        if Q > 1.0: Q = 1.0

                        levels = [2,4]
                        if TrophicComplexityLevel in levels and Rtype != 'c':
                            # b's are by-products of consuming a's
                            # c's are by-products of consuming b's

                            if Rtype == 'a':
                                R = R.replace('a', 'b')
                            if Rtype == 'b':
                                R = R.replace('b', 'c')

                            RList[j] = R
                            RVals[j] = R.count(Rtype)
                            RIDs[j] = RID
                            RID += 1

                        else:
                            # no resources produced as by-products of consuming other resources
                            # print 'reducing R'
                            RVals.pop(j)
                            RList.pop(j)
                            RIDs.pop(j)
                            RXs.pop(j)
                            RYs.pop(j)


                        Qs[i] = Q



    listlen = [len(SpeciesIDs), len(Qs), len(IndIDs), len(IXs), len(IYs), len(GrowthList), len(MaintList), len(DispList), len(ADList), len(EVList), len(TLList)]
    if min(listlen) != max(listlen):
        print 'In consume (bottom)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()

    return [RList, RVals, RIDs, RID, RXs, RYs, SpeciesIDs, Qs, encounters]





def reproduce(spec, SpeciesIDs, Qs, IDs, ID, Xs, Ys, w, h, GD, DispD, RD, MD, MFD,
    RPD, EnvD, envGs, GList, MList, DList, ADList, EVList, TLList, RList, RVals, RX, RY, RID, RIDs, TrophicComplexityLevel, SpatialComplexityLevel, ResourceComplexityLevel, BiologicalComplexityLevel):

    listlen = [len(SpeciesIDs), len(Qs), len(IDs), len(Xs), len(Ys), len(GList), len(MList), len(DList), len(ADList), len(EVList), len(TLList)]
    if min(listlen) != max(listlen):
        print 'In reproduce (top)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()


    if SpeciesIDs == []:
        return [SpeciesIDs, Qs, IDs, ID, Xs, Ys, GD, DispD, GList, MList, DList, ADList, EVList, TLList, RList, RVals, RX, RY, RID, RIDs]


    n = len(IDs)
    for j in range(n):

        if len(TLList) != len(EVList):
            print j, n, "len(TLList) != len(EVList) 2"
            print len(TLList), len(EVList)
            sys.exit()

        i = randint(0, len(IDs)-1)

        state = ADList[i]
        if state == 'd':
            continue

        Q = Qs[i]
        pq = float(np.mean(Q))

        if Q < 0.0:

            if TrophicComplexityLevel >= 3:
                RList.append('d')
                RVals.append(1)
                RX.append(Xs[i])
                RY.append(Ys[i])
                RIDs.append(RID)
                RID += 1

            Qs.pop(i)
            SpeciesIDs.pop(i)
            IDs.pop(i)
            Xs.pop(i)
            Ys.pop(i)
            GList.pop(i)
            MList.pop(i)
            TLList.pop(i)
            DList.pop(i)
            ADList.pop(i)
            EVList.pop(i)
            continue

        p = np.random.binomial(1, pq)
        if p == 1: # individual is large enough to reproduce

            spID = SpeciesIDs[i]
            X = Xs[i]
            Y = Ys[i]

            '''
            pg = []
            Speciesopts = EnvD[spID]

            for g, opt in enumerate(Speciesopts):

                x, y = envGs[g]
                pg.append(1 - (abs(X - x)/max([X, x])))
                pg.append(1 - (abs(Y - y)/max([Y, y])))


            if np.mean(pg) > 1 or np.mean(pg) < 0:
                print 'pg:', pg
                sys.exit()

            p = np.mean(pg)
            p = np.random.binomial(1, p)
            '''

            p = 1
            if p == 1: # the environment is suitable for reproduction

                Qs[i] = Q/2
                Qs.append(Q/2)

                tll = TLList[i]
                TLList.append(tll)

                ID += 1
                IDs.append(ID)

                p = np.random.binomial(1, spec)

                p = 0
                if p == 1:

                    # speciate
                    spID_new = float(spID)**3

                    # new species growth rate
                    p = np.random.binomial(1, 0.25)
                    if p == 1:
                        GD[spID_new] = np.random.uniform(0.5, 1.0)
                    else:
                        GD[spID_new] = GD[spID]

                    # new species maintenance
                    p = np.random.binomial(1, 0.25)
                    if p == 1:
                        MD[spID_new] = np.random.uniform(0.01, 0.1)
                    else:
                        MD[spID_new] = MD[spID]

                    # species environmental gradient optima
                    glist = []
                    for jj, g in enumerate(envGs):
                        p = np.random.binomial(1, 0.25)
                        if p == 1:
                            x = np.random.uniform(0.0, w)
                            y = np.random.uniform(0.0, h)
                        else:
                            x = EnvD[spID][jj][0]
                            y = EnvD[spID][jj][1]

                        glist.append([x, y])
                        EnvD[spID_new] = glist

                    # new species active dispersal rate
                    p = np.random.binomial(1, 0.25)
                    if p == 1: DispD[spID_new] = np.random.uniform(0.0, 0.1)
                    else: DispD[spID_new] = DispD[spID]

                    # new species resource use efficiency
                    p = np.random.binomial(1, 0.25)
                    if p == 1: RD[spID_new] = np.random.uniform(0.01, 0.9)
                    else: RD[spID_new] = RD[spID]

                    spID = spID_new

                means = GD[spID]
                ii = GetIndParam(means)
                GList.append(ii)

                means = MD[spID]
                ii = GetIndParam(means)
                MList.append(ii)

                means = RD[spID]
                ii = GetIndParam(means)
                EVList.append(ii)

                means = DispD[spID]
                ii = GetIndParam(means)
                DList.append(ii)

                SpeciesIDs.append(spID)
                ADList.append('a')

                newX = float(np.random.uniform(X-0.5, X+0.5, 1))
                if limit > newX: newX = 0
                if newX > w - limit: newX = w - limit
                Xs.append(newX)

                newY = float(np.random.uniform(Y-0.5, Y+0.5, 1))
                if limit > newY: newY = 0
                elif newY > h: newY = h - limit
                Ys.append(newY)


    listlen = [len(SpeciesIDs), len(Qs), len(IDs), len(Xs), len(Ys), len(GList), len(MList), len(DList), len(ADList), len(EVList), len(TLList)]
    if min(listlen) != max(listlen):
        print 'In reproduce (bottom)'
        print 'min(listlen) != max(listlen)'
        sys.exit()

    return [SpeciesIDs, Qs, IDs, ID, Xs, Ys, GD, DispD, GList, MList, DList, ADList, EVList, TLList, RList, RVals, RX, RY, RID, RIDs]






def res_dispersal(ct, RList, Vals, Xs, Ys, ID, IDs, numr, w, h, u0, TrophicComplexityLevel, SpatialComplexityLevel, ResourceComplexityLevel, BiologicalComplexityLevel):

    if IDs == []:
        return [RList, Vals, Xs, Ys, ID, IDs]

    n = len(IDs)
    for j in range(n):

        i = randint(0, len(IDs)-1)

        Ys[i] = float(np.random.uniform(0.001*h, 0.999*h))
        Xs[i] = float(np.random.uniform(0.001*w, 0.999*w))

    return [RList, Vals, Xs, Ys, ID, IDs]





def dispersal(spec, SpeciesIDs, Qs, IDs, ID, Xs, Ys,  w, h, GD, DispD, RD, MD, MFD,
        RPD, EnvD, envGs, GList, MList, DList, ADList, EVList, TLList, RList, RVals, RXs, RYs, RIDs, RID, TrophicComplexityLevel, SpatialComplexityLevel, ResourceComplexityLevel, BiologicalComplexityLevel):

    listlen = [len(SpeciesIDs), len(Qs), len(IDs), len(Xs), len(Ys), len(GList), len(MList), len(DList), len(ADList), len(EVList), len(TLList)]
    if min(listlen) != max(listlen):
        print 'In dispersal (top)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()

    if SpeciesIDs == []:
        return [SpeciesIDs, Qs, IDs, ID, Xs, Ys, GD, DispD, GList, MList, DList, ADList, EVList, TLList, RList, RVals, RXs, RYs, RIDs, RID]

    n = len(IDs)
    for j in range(n):

        i = randint(0, len(IDs)-1)

        spID = SpeciesIDs[i]
        X = Xs[i]
        Y = Ys[i]


        if SpatialComplexityLevel == 1:
            Ys[i] = float(np.random.uniform(0.001*h, 0.999*h))
            Xs[i] = float(np.random.uniform(0.001*w, 0.999*w))

        state = ADList[i]
        if state == 'd':
            continue

        if SpatialComplexityLevel == 2:
            dist = DispD[spID]

            r = Qs[i]

            if r < 0.0:
                print 'line 653: Q:', r
                sys.exit()

            if Qs[i] >= MList[i]*dist:

                # A cost for active dispersal
                if BiologicalComplexityLevel > 1:
                    r -= MList[i]*dist

                    if r < 0.0:

                        if TrophicComplexityLevel >= 3:
                            RList.append('d')
                            RVals.append(1)
                            RXs.append(Xs[i])
                            RYs.append(Ys[i])
                            RIDs.append(RID)
                            RID += 1

                        Qs.pop(i)
                        SpeciesIDs.pop(i)
                        IDs.pop(i)
                        Xs.pop(i)
                        Ys.pop(i)
                        GList.pop(i)
                        MList.pop(i)
                        TLList.pop(i)
                        DList.pop(i)
                        ADList.pop(i)
                        EVList.pop(i)
                        continue

                    Qs[i] = r

                vd = choice([-1, 1])
                hd = choice([-1, 1])

                X += hd*dist
                Y += vd*dist

                if X > 0.999 * w: X = w * 0.999
                elif X < 0.001 * w: X = 0.001 * w
                if Y > 0.999 * h: Y = h * 0.999
                elif Y < 0.001 * h: Y = 0.001 * h

                Xs[i] = X
                Ys[i] = Y

    listlen = [len(SpeciesIDs), len(Qs), len(IDs), len(Xs), len(Ys), len(GList), len(MList), len(DList), len(ADList), len(EVList), len(TLList)]
    if min(listlen) != max(listlen):
        print 'In dispersal (bottom)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()


    return [SpeciesIDs, Qs, IDs, ID, Xs, Ys, GD, DispD, GList, MList, DList, ADList, EVList, TLList, RList, RVals, RXs, RYs, RIDs, RID]


def chemotaxis(RList, RVals, RIDs, RID, RXs, RYs, SpeciesIDs, Qs, IndIDs, IndID,
    IXs, IYs, w, h, GD, RD, DispD, GrowthList, MaintList, DispList, ADList, TLList, EVList,
    TrophicComplexityLevel, SpatialComplexityLevel, ResourceComplexityLevel, BiologicalComplexityLevel, enz_action = False):

    listlen = [len(SpeciesIDs), len(Qs), len(IndIDs), len(IXs), len(IYs), len(GrowthList), len(MaintList), len(DispList), len(ADList), len(EVList), len(TLList)]

    if min(listlen) != max(listlen):
        print 'In chemotaxis (top)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()

    if SpeciesIDs == [] or RIDs == []:
        return [RList, RVals, RIDs, RID, RXs, RYs, SpeciesIDs, Qs, IndIDs, IndID, IXs, IYs, w, h, GD, RD, DispD, GrowthList, MaintList, DispList, ADList, TLList, EVList]

    n = len(IndIDs)

    for ii in range(n):

        n = len(IndIDs)
        i = randint(0, n-1)

        state = ADList[i]
        if state == 'd':
            continue

        spID = SpeciesIDs[i]
        disp = DispD[spID]

        tl = TLList[i] # Trophic level
        x1 = IXs[i]
        y1 = IYs[i]

        Try = min([50, len(RIDs)])
        ct = 0

        minDist = 100
        targetX = 0
        targetY = 0
        dist = 0

        while ct < Try and RIDs is not []:

            ct += 1

            r = len(RIDs)
            Try = min([50, r])
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

                if dist < minDist:
                    minDist = dist
                    targetX = RXs[j]
                    targetY = RYs[j]

        if BiologicalComplexityLevel > 1:

            # A cost for active dispersal
            r = Qs[i]
            if r - MaintList[i]*disp > MaintList[i]:
                r -= MaintList[i]*disp
            else:
                dist = (r - MaintList[i])/MaintList[i]
                r -= MaintList[i]*disp

            if r < 0.0:

                if TrophicComplexityLevel >= 3:
                    RList.append('d')
                    RVals.append(1)
                    RXs.append(IXs[i])
                    RYs.append(IYs[i])
                    RIDs.append(RID)
                    RID += 1

                Qs.pop(i)
                SpeciesIDs.pop(i)
                IndIDs.pop(i)
                IXs.pop(i)
                IYs.pop(i)
                GrowthList.pop(i)
                MaintList.pop(i)
                TLList.pop(i)
                DispList.pop(i)
                ADList.pop(i)
                EVList.pop(i)
                continue

            Qs[i] = r

        if x1 > targetX:
            x1 -= dist*disp

        elif x1 < targetX:
            x1 += dist*disp

        if y1 > targetY:
            y1 -= dist*disp

        elif y1 < targetY:
            y1 += dist*disp

        if x1 > 0.999 * w: x1 = w * 0.999
        elif x1 < 0.001 * w: x1 = 0.001 * w
        if y1 > 0.999 * h: y1 = h * 0.999
        elif y1 < 0.001 * h: y1 = 0.001 * h

        IXs[i] = x1
        IYs[i] = y1


    listlen = [len(SpeciesIDs), len(Qs), len(IndIDs), len(IXs), len(IYs), len(GrowthList), len(MaintList), len(DispList), len(ADList), len(EVList), len(TLList)]
    if min(listlen) != max(listlen):
        print 'In chemotaxis (bottom)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()

    return [RList, RVals, RIDs, RID, RXs, RYs, SpeciesIDs, Qs, IndIDs, IndID, IXs, IYs, w, h, GD, RD, DispD, GrowthList, MaintList, DispList, ADList, TLList, EVList]
