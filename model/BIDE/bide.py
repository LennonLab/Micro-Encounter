# -*- coding: utf-8 -*-
from __future__ import division
from random import randint, choice
import numpy as np
import sys
import math
#from math import modf
#import decimal
import time


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
        RAD.append(vector.count(val)) # the abundance of each Sp_

    return RAD, unique # the rad and the specieslist


def ResIn(RList, Vals, Xs, Ys, ID, IDs, numr, rmax, w, h, u0):

    for r in range(numr):
        x = np.random.binomial(1, u0)

        if x == 1:
            rval = int(np.random.random_integers(1, rmax, 1))
            Vals.append(rval)

            IDs.append(ID)
            ID += 1

            r = choice(['A', 'B', 'C'])
            RList.append(r)

            Ys.append(float(np.random.uniform(0.1*h, 0.9*h)))
            Xs.append(float(np.random.uniform(0.1*w, 0.9*w)))

    return [RList, Vals, Xs, Ys, IDs, ID]



def immigration(mfmax, p_max, d_max, g_max, m_max, seed, ip, Sp, Xs, Ys, w, h, MD, MFD, RPD,
        EnvD, envGs, GD, DispD, IDs, ID, Qs, RD, u0, alpha, GList, MList, DList, ADList, EVList, TLList, ct):

    listlen = [len(Sp), len(Qs), len(IDs), len(Xs), len(Ys), len(GList), len(MList), len(DList), len(ADList)]
    if min(listlen) != max(listlen):
        print 'In immigration (top)'
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

            Ys.append(float(np.random.uniform(0.1*h, 0.9*h)))
            Xs.append(float(np.random.uniform(0.1*w, 0.9*w)))

            IDs.append(ID)
            ID += 1

            tl = choice(['A', 'B', 'C'])
            TLList.append(tl)

            Q = float(np.random.uniform(0.05, 0.5))
            Qs.append(Q)

            if prop not in GD:

                # species growth rate
                GD[prop] = np.random.uniform(g_max/10, g_max)

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
                for g in envGs:
                    x = np.random.uniform(0.0, w)
                    y = np.random.uniform(0.0, h)
                    glist.append([x,y])
                EnvD[prop] = glist

                # Resource use efficiency
                RD[prop] = np.random.uniform(0.1, 1.0, 1)

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


            listlen = [len(Sp), len(Qs), len(IDs), len(Xs), len(Ys), len(GList), len(MList), len(DList), len(ADList)]
            if min(listlen) != max(listlen):
                print 'In immigration'
                print 'min(listlen) != max(listlen)'
                print listlen
                sys.exit()


    return [Sp, Xs, Ys, MD, MFD, RPD, EnvD, GD, DispD, IDs, ID, Qs,
            RD, GList, MList, DList, ADList, EVList, TLList]





def maintenance(Sp_IDs, Xs, Ys, MD, MFD, RPD, EnvD, IDs, Qs, GrowthList,
        MaintList, TLList, DList, ADList):

    if Sp_IDs == []:
        return [Sp_IDs, Xs, Ys, IDs, Qs, GrowthList, MaintList, TLList, DList, ADList]

    n = len(IDs)
    for j in range(n):

        i = randint(0, len(IDs)-1)

        val = Qs[i]
        val -= MaintList[i] # maintanence influenced by species id

        if val <= MaintList[i]*0.00001:   # starved

            Qs.pop(i)
            Sp_IDs.pop(i)
            IDs.pop(i)
            Xs.pop(i)
            Ys.pop(i)
            GrowthList.pop(i)
            MaintList.pop(i)
            TLList.pop(i)
            DList.pop(i)
            ADList.pop(i)

        else: Qs[i] = val

    return [Sp_IDs, Xs, Ys, IDs, Qs, GrowthList, MaintList, TLList, DList, ADList]





def transition(Sp_IDs, IDs, Qs, GrowthList, MaintList, MFD, RPD, ADList):

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

    return [Sp_IDs, IDs, Qs, GrowthList, MaintList, ADList]





def consume(field, RList, R_Vals, R_IDs, R_ID, RXs, RYs, Sp_IDs, Qs, I_IDs, I_ID,
    IXs, IYs, w, h, GD, RD, DispD, GrowthList, MaintList, DispList, ADList, TLList, EVList, enz_action = True):

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

            #print tl, R

            if tl == R:
                r = str()

                if R == 'A':
                    R = 'B'
                    r = 'B'
                elif R == 'B':
                    R = 'C'
                    r = 'C'

                x2 = RXs[j]
                y2 = RYs[j]
                dist = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

                ind_radius = np.mean(Qs[i])
                R_val = R_Vals[j]
                res_radius = R_val

                if dist <= ind_radius + res_radius:

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

    return [RList, R_Vals, R_IDs, R_ID, RXs, RYs, Sp_IDs, Qs]





def reproduce(spec, Sp_IDs, Qs, IDs, ID, Xs, Ys, w, h, GD, DispD, RD, MD, MFD,
    RPD, EnvD, envGs, GList, MList, RList, DList, ADList):

    listlen = [len(Sp_IDs), len(Qs), len(IDs), len(Xs), len(Ys), len(GList), len(MList), len(DList), len(ADList)]
    if min(listlen) != max(listlen):
        print 'In reproduce'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()


    if Sp_IDs == []:
        return [Sp_IDs, Qs, IDs, ID, Xs, Ys, GD, DispD, GList, MList, RList, DList, ADList]

    n = len(IDs)
    for j in range(n):

        i = randint(0, len(IDs)-1)

        state = ADList[i]
        if state == 'd':
            continue

        Q = Qs[i]
        pq = float(np.mean(Q))
        p = np.random.binomial(1, pq)

        if p == 1: # individual is large enough to reproduce

            spID = Sp_IDs[i]
            X = Xs[i]
            Y = Ys[i]

            # A cost for reproducing
            r = Qs[i]
            r -= MList[i]*r
            Qs[i] = r

            pg = []
            sp_opts = EnvD[spID]

            for g, opt in enumerate(sp_opts):

                x, y = envGs[g]
                pg.append(1 - (abs(X - x)/max([X, x])))
                pg.append(1 - (abs(Y - y)/max([Y, y])))


            if np.mean(pg) > 1 or np.mean(pg) < 0:
                print pg
                sys.exit()

            p = np.mean(pg)
            p = np.random.binomial(1, p)
            if p == 1: # the environment is suitable for reproduction

                Qs[i] = Q/2
                Qs.append(Q/2)

                ID += 1
                IDs.append(ID)

                p = np.random.binomial(1, spec)
                p = 0
                if p == 1:

                    # speciate
                    spID_new = spID**3

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
                    for j, g in enumerate(envGs):
                        p = np.random.binomial(1, 0.25)
                        if p == 1:
                            x = np.random.uniform(0.0, w)
                            y = np.random.uniform(0.0, h)
                        else:
                            x = EnvD[spID][j][0]
                            y = EnvD[spID][j][1]

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
                i = GetIndParam(means)
                GList.append(i)

                means = MD[spID]
                i = GetIndParam(means)
                MList.append(i)

                means = RD[spID]
                i = GetIndParam(means)
                RList.append(i)

                means = DispD[spID]
                i = GetIndParam(means)
                DList.append(i)

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

    return [Sp_IDs, Qs, IDs, ID, Xs, Ys, GD, DispD, GList, MList, RList, DList, ADList]




def chemotaxis(spec, Sp_IDs, Qs, IDs, ID, Xs, Ys,  w, h, GD, DispD, RD, MD, MFD,
        RPD, EnvD, envGs, GList, MList, RList, DList, ADList):

    if Sp_IDs == []:
        return [Sp_IDs, Qs, IDs, ID, Xs, Ys, GD, DispD, GList, MList, RList, DList, ADList]

    n = len(IDs)
    for j in range(n):

        i = randint(0, len(IDs)-1)

        state = ADList[i]
        if state == 'd':
            continue

        spID = Sp_IDs[i]
        X = Xs[i]
        Y = Ys[i]

        sp_opts = EnvD[spID]

        for g, opt in enumerate(sp_opts):
            x, y = opt
            dist = DispD[spID]

            if g == 0:
                # A cost for active dispersal

                r = Qs[i]
                r -= MList[i]*dist*r
                Qs[i] = r

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


    return [Sp_IDs, Qs, IDs, ID, Xs, Ys, GD, DispD, GList, MList, RList, DList, ADList]
