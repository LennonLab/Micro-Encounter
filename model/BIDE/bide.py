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


def randSplit(word, splits):
    for splitLen in splits:
        if splitLen > len(word):
            break
        yield word[:splitLen]
        word = word[splitLen::]


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
        RAD.append(vector.count(val)) # the abundance of each Sp_

    return RAD, unique # the rad and the specieslist


def ResIn(ct, RList, Vals, Xs, Ys, ID, IDs, numr, rmax, w, h, u0, TrophicComplexityLevel, SpatialComplexityLevel, ResourceComplexityLevel, BiologicalComplexityLevel):

    Ymean = float(np.random.uniform(0.1*h, 0.9*h))
    Xmean = float(np.random.uniform(0.1*w, 0.9*w))

    listlen = [len(RList), len(Vals), len(Xs), len(Ys), len(IDs)]
    if min(listlen) != max(listlen):
        print ct, 'In ResIn (top)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()


    for r in range(numr):
        x = np.random.binomial(1, u0)

        if x == 1:
            rval = int(np.random.random_integers(1, rmax, 1))
            Vals.append(rval)

            IDs.append(ID)
            ID += 1

            if ResourceComplexityLevel == 1:
                r = 'a'

            elif ResourceComplexityLevel == 2:
                r = choice(['a', 'b', 'c'])

            elif ResourceComplexityLevel == 3:

                a = ['a-a-a-a-a-a-a-a-a-a-a-a-a-a-a-a-a-a-a-a-a-a-a-a-a-a-a-a-a-a-a-a']
                b = ['bbb-bbb-bbb-bbb-bbb-bbb-bbb-bbb-bbb-bbb-bbb-bbb-bbb-bbb-bbb-bbb']
                c = ['cccccc-cccccc-cccccc-cccccc-cccccc-cccccc-cccccc-cccccc-cccccc']
                r = choice([a, b, c])

            RList.append(r)

            if SpatialComplexityLevel < 3:

                Ys.append(float(np.random.uniform(0.1*h, 0.9*h)))
                Xs.append(float(np.random.uniform(0.1*w, 0.9*w)))

            elif SpatialComplexityLevel == 3:

                std = 0.01
                vals = np.random.normal([Ymean, Xmean], std)

                if vals[0] < 0.1*h:
                    vals[0] = 0.1*h
                elif vals[0] > 0.9*h:
                    vals[0] = 0.9*h

                if vals[1] < 0.1*w:
                    vals[1] = 0.1*w
                elif vals[1] > 0.9*w:
                    vals[1] = 0.9*w

                Ys.append(vals[0])
                Xs.append(vals[1])

    listlen = [len(RList), len(Vals), len(Xs), len(Ys), len(IDs)]
    if min(listlen) != max(listlen):
        print 'In ResIn (bottom)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()

    return [RList, Vals, Xs, Ys, IDs, ID]



def immigration(mfmax, p_max, d_max, g_max, m_max, seed, ip, Sp, Xs, Ys, w, h, MD,
    MFD, RPD, EnvD, envGs, GD, DispD, IDs, ID, Qs, RD, u0, alpha, GList, MList,
    DList, ADList, EVList, TLList, ct, TrophicComplexityLevel, SpatialComplexityLevel,
    ResourceComplexityLevel, BiologicalComplexityLevel):

    Ymean = float(np.random.uniform(0.1*h, 0.9*h))
    Xmean = float(np.random.uniform(0.1*w, 0.9*w))

    listlen = [len(Sp), len(Qs), len(IDs), len(Xs), len(Ys), len(GList), len(MList), len(DList), len(ADList), len(EVList), len(TLList)]
    if min(listlen) != max(listlen):
        print ct, 'In immigration (top)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()


    if ct > 1: seed = 1

    for m in range(seed):
        x = 0

        if seed > 1: x = 1

        else: x = np.random.binomial(1, u0*ip)

        if x == 1:
            prop = str(float(np.random.logseries(alpha, 1)))

            Sp.append(prop)

            if SpatialComplexityLevel < 3:
                Ys.append(float(np.random.uniform(0.1*h, 0.9*h)))
                Xs.append(float(np.random.uniform(0.1*w, 0.9*w)))

            elif SpatialComplexityLevel == 3:

                std = 0.01
                vals = np.random.normal([Ymean, Xmean], std)

                if vals[0] < 0.1*h:
                    vals[0] = 0.1*h
                elif vals[0] > 0.9*h:
                    vals[0] = 0.9*h

                if vals[1] < 0.1*w:
                    vals[1] = 0.1*w
                elif vals[1] > 0.9*w:
                    vals[1] = 0.9*w

                Ys.append(vals[0])
                Xs.append(vals[1])

            IDs.append(ID)
            ID += 1

            tl = choice(['a', 'b', 'c'])
            TLList.append(tl)

            Q = float(np.random.uniform(0.05, 0.5))
            Qs.append(Q)

            if prop not in GD:

                # species growth rate
                g = np.random.uniform(g_max/10, g_max)
                GD[prop] = g

                # species maintenance
                MD[prop] = np.random.uniform(m_max/10, m_max)

                # species maintenance reduction factor
                MFD[prop] = np.random.uniform(1, mfmax)

                # species resuscitation factor
                RPD[prop] = np.random.uniform(0.001, 0.1)

                # species active dispersal rate
                DispD[prop] = np.random.uniform(d_max/10, d_max)

                # species environmental gradient optima
                glist = []
                for e in envGs:
                    x = np.random.uniform(0.0, w)
                    y = np.random.uniform(0.0, h)
                    glist.append([x,y])
                EnvD[prop] = glist

                # Resource use efficiency
                RD[prop] = [1-g]

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






def maintenance(Sp_IDs, Xs, Ys, MD, MFD, RPD, EnvD, IDs, Qs, GList, MaintList,
    DList, ADList, EVList, TLList, TrophicComplexityLevel, SpatialComplexityLevel,
    ResourceComplexityLevel, BiologicalComplexityLevel):

    listlen = [len(Sp_IDs), len(Qs), len(IDs), len(Xs), len(Ys), len(GList), len(MaintList), len(DList), len(ADList), len(EVList), len(TLList)]
    if min(listlen) != max(listlen):
        print 'In maintenance (top)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()

    if Sp_IDs == []:
        return [Sp_IDs, Xs, Ys, IDs, Qs, GList, MaintList, DList, ADList, EVList, TLList]

    n = len(IDs)
    for j in range(n):

        i = randint(0, len(IDs)-1)

        val = Qs[i]
        val -= MaintList[i] # maintanence influenced by species id

        if val <= MaintList[i]*0.001:   # starved

            Qs.pop(i)
            Sp_IDs.pop(i)
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

    listlen = [len(Sp_IDs), len(Qs), len(IDs), len(Xs), len(Ys), len(GList), len(MaintList), len(DList), len(ADList), len(EVList), len(TLList)]
    if min(listlen) != max(listlen):
        print 'In maintenance (bottom)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()

    return [Sp_IDs, Xs, Ys, IDs, Qs, GList, MaintList, DList, ADList, EVList, TLList]





def transition(Sp_IDs, IDs, Qs, GrowthList, MaintList, MFD, RPD, ADList, TrophicComplexityLevel, SpatialComplexityLevel, ResourceComplexityLevel, BiologicalComplexityLevel):

    if Sp_IDs == []:
        return [Sp_IDs, IDs, Qs, GrowthList, MaintList, ADList]

    n = len(IDs)
    for j in range(n):

        i = randint(0, len(IDs)-1)
        spid = Sp_IDs[i]
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
                print 'line 321: Q:', Q
                sys.exit()


    return [Sp_IDs, IDs, Qs, GrowthList, MaintList, ADList]





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
            R = RList[j][0]
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
                    ct = Try

                    if '-' in R:
                        # useful stackoverflow answers:
                        # http://stackoverflow.com/questions/32212998/split-string-by-n-when-n-is-random
                        # http://stackoverflow.com/questions/19954593/python-checking-a-strings-first-and-last-character
                        # http://stackoverflow.com/questions/5188792/how-to-check-a-string-for-specific-characters

                        pos = randint(0, len(R))
                        eat = 'no'

                        if R[pos] == '-':
                            pieces = randSplit(R, randint(0,len(R)))

                            for index, piece in enumerate(pieces):
                                if piece.startswith('-'):
                                    piece = piece[1:]

                                if piece.endswith('-'):
                                    piece = piece[:-1]

                            if '-' in piece[0] and '-' in piece[1]:
                                RList[j] = pieces[0]
                                R_Vals[j] = pieces[0].count(Rtype)

                                RList.append(pieces[1])
                                R_Vals[j] = pieces[1].count(Rtype)
                                R_IDs.append(R_ID)
                                R_ID += 1
                                RXs.append(x2)
                                RYs.append(y2)


                            elif '-' not in piece[0]:
                                eat = 'yes'
                                RList[j] = pieces[1]
                                R_Vals[j] = pieces[1].count(Rtype)

                                R = piece[0]
                                R_val = pieces[0].count(Rtype)

                            elif '-' not in piece[1]:
                                eat = 'yes'
                                RList[j] = pieces[0]
                                R_Vals[j] = pieces[0].count(Rtype)

                                R = piece[1]
                                R_val = pieces[1].count(Rtype)

                            if eat == 'yes':

                                # The individual's cell quota
                                Q = Qs[i]

                                # The species
                                sp = Sp_IDs[i]
                                mu = GD[sp]
                                efficiency = EVList[i][0]
                                mu = mu * efficiency

                                # Increase cell quota & decrease resource particle size
                                if R_val > (mu * Q):
                                    R_val = R_val - (mu * Q)*10
                                    Q += (mu * Q)

                                else:
                                    Q += R_val
                                    R_val = 0.0

                                if Q > 1.0:
                                    R_val = R_val + (Q - 1.0)
                                    Q = 1.0
                                    R_Vals[j] = R_val

                                if R_val <= 0.0:
                                    R_Vals.pop(j)
                                    RList.pop(j)
                                    R_IDs.pop(j)
                                    RXs.pop(j)
                                    RYs.pop(j)

                                Qs[i] = Q


    listlen = [len(Sp_IDs), len(Qs), len(I_IDs), len(IXs), len(IYs), len(GrowthList), len(MaintList), len(DispList), len(ADList), len(EVList), len(TLList)]
    if min(listlen) != max(listlen):
        print 'In consume (bottom)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()

    return [RList, R_Vals, R_IDs, R_ID, RXs, RYs, Sp_IDs, Qs]





def reproduce(spec, Sp_IDs, Qs, IDs, ID, Xs, Ys, w, h, GD, DispD, RD, MD, MFD,
    RPD, EnvD, envGs, GList, MList, DList, ADList, EVList, TLList, TrophicComplexityLevel, SpatialComplexityLevel, ResourceComplexityLevel, BiologicalComplexityLevel):

    listlen = [len(Sp_IDs), len(Qs), len(IDs), len(Xs), len(Ys), len(GList), len(MList), len(DList), len(ADList), len(EVList), len(TLList)]
    if min(listlen) != max(listlen):
        print 'In reproduce (top)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()


    if Sp_IDs == []:
        return [Sp_IDs, Qs, IDs, ID, Xs, Ys, GD, DispD, GList, MList, DList, ADList, EVList, TLList]


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
            print 'line 456: Q:', Q
            sys.exit()

        p = np.random.binomial(1, pq)

        p = 1
        if p == 1: # individual is large enough to reproduce

            spID = Sp_IDs[i]
            X = Xs[i]
            Y = Ys[i]

            pg = []
            sp_opts = EnvD[spID]

            for g, opt in enumerate(sp_opts):

                x, y = envGs[g]
                pg.append(1 - (abs(X - x)/max([X, x])))
                pg.append(1 - (abs(Y - y)/max([Y, y])))


            if np.mean(pg) > 1 or np.mean(pg) < 0:
                print 'pg:', pg
                sys.exit()

            p = np.mean(pg)
            p = np.random.binomial(1, p)

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

                Sp_IDs.append(spID)
                ADList.append('a')

                newX = float(np.random.uniform(X-0.5, X+0.5, 1))
                if limit > newX: newX = 0
                if newX > w - limit: newX = w - limit
                Xs.append(newX)

                newY = float(np.random.uniform(Y-0.5, Y+0.5, 1))
                if limit > newY: newY = 0
                elif newY > h: newY = h - limit
                Ys.append(newY)


    listlen = [len(Sp_IDs), len(Qs), len(IDs), len(Xs), len(Ys), len(GList), len(MList), len(DList), len(ADList), len(EVList), len(TLList)]
    if min(listlen) != max(listlen):
        print 'In reproduce (bottom)'
        print 'min(listlen) != max(listlen)'
        sys.exit()

    return [Sp_IDs, Qs, IDs, ID, Xs, Ys, GD, DispD, GList, MList, DList, ADList, EVList, TLList]






def res_dispersal(ct, RList, Vals, Xs, Ys, ID, IDs, numr, rmax, w, h, u0, TrophicComplexityLevel, SpatialComplexityLevel, ResourceComplexityLevel, BiologicalComplexityLevel):

    if IDs == []:
        return [RList, Vals, Xs, Ys, ID, IDs]

    n = len(IDs)
    for j in range(n):

        i = randint(0, len(IDs)-1)

        Ys[i] = float(np.random.uniform(0.1*h, 0.9*h))
        Xs[i] = float(np.random.uniform(0.1*w, 0.9*w))

    return [RList, Vals, Xs, Ys, ID, IDs]





def dispersal(spec, Sp_IDs, Qs, IDs, ID, Xs, Ys,  w, h, GD, DispD, RD, MD, MFD,
        RPD, EnvD, envGs, GList, MList, DList, ADList, EVList, TLList, TrophicComplexityLevel, SpatialComplexityLevel, ResourceComplexityLevel, BiologicalComplexityLevel):

    listlen = [len(Sp_IDs), len(Qs), len(IDs), len(Xs), len(Ys), len(GList), len(MList), len(DList), len(ADList), len(EVList), len(TLList)]
    if min(listlen) != max(listlen):
        print 'In dispersal (top)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()

    if Sp_IDs == []:
        return [Sp_IDs, Qs, IDs, ID, Xs, Ys, GD, DispD, GList, MList, DList, ADList, EVList, TLList]

    n = len(IDs)
    for j in range(n):

        i = randint(0, len(IDs)-1)

        spID = Sp_IDs[i]
        X = Xs[i]
        Y = Ys[i]


        if SpatialComplexityLevel == 1:
            Ys[i] = float(np.random.uniform(0.1*h, 0.9*h))
            Xs[i] = float(np.random.uniform(0.1*w, 0.9*w))

        state = ADList[i]
        if state == 'd':
            continue

        if SpatialComplexityLevel == 2:
            dist = DispD[spID]*10

            r = Qs[i]

            if r < 0.0:
                print 'line 653: Q:', r
                sys.exit()

            if Qs[i] >= MList[i]*dist:

                # A cost for active dispersal
                r -= MList[i]*dist
                Qs[i] = r

                vd = choice([-1, 1])
                hd = choice([-1, 1])

                X += hd*dist
                Y += vd*dist

                if X > w: X = w
                elif X < 0: X = 0
                if Y > h: Y = h
                elif Y < 0: Y = 0

                Xs[i] = X
                Ys[i] = Y

        elif SpatialComplexityLevel == 3:
            sp_opts = EnvD[spID]

            move = 'n'
            for g, opt in enumerate(sp_opts):
                x, y = opt
                dist = DispD[spID]

                if g == 0:

                    # A cost for active dispersal
                    if Qs[i] >= MList[i]*dist:

                        r = Qs[i]
                        r -= MList[i]*dist
                        if r < 0.0:
                            print 'line 689: Q:', r
                            sys.exit()

                        Qs[i] = r
                        move = 'y'

                if move == 'y':
                    if x > X:
                        X += dist

                    elif x < X:
                        X -= dist

                    if y > Y:
                        Y += dist

                    elif y < Y:
                        Y -= dist

                    if X > w: X = w
                    elif X < 0: X = 0
                    if Y > h: Y = h
                    elif Y < 0: Y = 0

                Xs[i] = X
                Ys[i] = Y

    listlen = [len(Sp_IDs), len(Qs), len(IDs), len(Xs), len(Ys), len(GList), len(MList), len(DList), len(ADList), len(EVList), len(TLList)]
    if min(listlen) != max(listlen):
        print 'In dispersal (bottom)'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()


    return [Sp_IDs, Qs, IDs, ID, Xs, Ys, GD, DispD, GList, MList, DList, ADList, EVList, TLList]
