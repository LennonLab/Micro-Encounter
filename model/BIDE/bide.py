# -*- coding: utf-8 -*-
from __future__ import division
from random import randint, choice
import numpy as np
import sys
import math
#from math import modf
#import decimal
import time

sys.path.append(mydir + "GitHub/Micro-Encounter/model/spatial")
import spatial


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


def ResIn(Type, Vals, Xs, Ys, ID, IDs, numr, rmax, nN, nP, nC, w, h, u0):


    for r in range(numr):
        x = np.random.binomial(1, u0)

        if x == 1:
            rval = int(np.random.random_integers(1, rmax, 1))
            nr = choice(['N', 'P', 'C'])

            if nr == 'N':
                rtype = int(np.random.random_integers(0, nN-1, 1))
                rtype = 'N'+str(rtype)

            if nr == 'P':
                rtype = int(np.random.random_integers(0, nP-1, 1))
                rtype = 'P'+str(rtype)

            if nr == 'C':
                rtype = int(np.random.random_integers(0, nC-1, 1))
                rtype = 'C'+str(rtype)

            Vals.append(rval)
            IDs.append(ID)
            Type.append(rtype)
            ID += 1

            Ys.append(float(np.random.uniform(0.1*h, 0.9*h)))
            Xs.append(float(np.random.uniform(0.1*w, 0.9*w)))

    return [Type, Vals, Xs, Ys, IDs, ID]



def immigration(mfmax, p_max, d_max, g_max, m_max, seed, ip, Sp, Xs, Ys, w, h, MD, MFD, RPD,
        EnvD, envGs, GD, DispD, IDs, ID, Qs, N_RD, P_RD, C_RD, nN, nP, nC, u0, alpha, GList,
        MList, NList, PList, CList, DList, ADList, ct):

    if u0 > 1.0:
        u0 = 1.0

    for m in range(seed):
        x = 0

        if seed > 1:
            x = 1

        else:
            x = np.random.binomial(1, u0*ip)

        if x == 1:
            prop = str(float(np.random.logseries(alpha, 1)))

            Sp.append(prop)

            Ys.append(float(np.random.uniform(0.1*h, 0.9*h)))
            Xs.append(float(np.random.uniform(0.1*w, 0.9*w)))

            IDs.append(ID)
            ID += 1

            Qn = float(np.random.uniform(0.01, 0.1))
            Qp = float(np.random.uniform(0.01, 0.1))
            Qc = float(np.random.uniform(0.01, 0.1))
            Qs.append([Qn, Qp, Qc])

            if prop not in GD:
                # species growth rate
                GD[prop] = np.random.uniform(g_max/10, g_max)

                # species maintenance
                MD[prop] = np.random.uniform(m_max/10, m_max)

                # species maintenance factor
                MFD[prop] = np.random.uniform(1, mfmax)

                # species RPF factor
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

                # species Nitrogen use efficiency
                N_RD[prop] = np.random.uniform(0.1, 1.0, nN)

                # species Phosphorus use efficiency
                P_RD[prop] = np.random.uniform(0.1, 1.0, nP)

                # species Carbon use efficiency
                C_RD[prop] = np.random.uniform(0.1, 1.0, nC)

            state = choice(['a','d'])
            ADList.append(state)

            means = GD[prop]
            i = GetIndParam(means)

            if state == 'a':
                GList.append(i)
            elif state == 'd':
                GList.append(i)

            means = MD[prop]
            i = GetIndParam(means)

            if state == 'a':
                MList.append(i)
            if state == 'd':
                MList.append(i/2.0)

            means = N_RD[prop]
            n = GetIndParam(means)
            means = P_RD[prop]
            p = GetIndParam(means)
            means = C_RD[prop]
            c = GetIndParam(means)

            NList.append(n)
            PList.append(p)
            CList.append(c)


            means = DispD[prop]
            i = GetIndParam(means)
            DList.append(i)


    return [Sp, Xs, Ys, MD, MFD, RPD, EnvD, GD, DispD, IDs, ID, Qs, N_RD,
            P_RD, C_RD, GList, MList, NList, PList, CList, DList, ADList]




def movement(TypeOf, motion, List, Xs, Ys, w, h, u0):

    Type, IDs, ID, Vals = [], [], int(), []

    if TypeOf == 'resource':
        Type, IDs, ID, Vals = List
    elif TypeOf == 'individual':
        Type, IDs, ID, Vals, DispD, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList, ADList = List

    else:
        IDs = List

    if Xs == []:
        if TypeOf == 'tracer':
            return [IDs, Xs, Ys]
        elif TypeOf == 'individual':
            return [Type, Xs, Ys, IDs, ID, Vals, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList, ADList]
        elif TypeOf == 'resource':
            return [Type, Xs, Ys, IDs, ID, Vals]


    limit, distance, direction, pop = 0.1, 0, 0, 'no'
    x, y = 0, 0

    n = len(IDs)
    for j in range(n):

        i = randint(0, len(IDs)-1)

        # get distance
        if TypeOf == 'individual':
            distance = np.random.uniform(0, DispList[i])
        else:
            distance = np.random.uniform(0, u0)

        x, y = Xs[i], Ys[i]

        # Go up or down
        if motion == 'unidirectional':
            direction = 1

        elif motion == 'brown_noise':
            direction = choice([-1, 1])

        y = y + (direction * distance)


        # get distance
        if TypeOf == 'individual':
            distance = np.random.uniform(0, DispList[i])
        else:
            distance = np.random.uniform(0, u0)

        # go forward or backward
        if TypeOf == 'resource':
            direction = 0.0
        else:
            if motion == 'unidirectional':
                direction = 1
            if motion == 'brown_noise':
                direction = choice([-1, 1])

        x = x + (direction * distance)

        if x > w - limit or x < limit:
            pop = 'yes'
        elif y > h - limit or y < limit:
            pop = 'yes'


        if motion == 'white_noise':
            y = float(np.random.uniform(0.1*h, 0.9*h))
            x = float(np.random.uniform(0.1*w, 0.9*w))

        if pop == 'no':
            Xs[i], Ys[i] = x, y

        elif pop == 'yes':
            Xs.pop(i)
            Ys.pop(i)
            IDs.pop(i)

            if TypeOf == 'resource' or TypeOf == 'individual':
                Type.pop(i)
                Vals.pop(i)

            if TypeOf == 'individual':
                GrowthList.pop(i)
                MaintList.pop(i)
                N_RList.pop(i)
                P_RList.pop(i)
                C_RList.pop(i)
                DispList.pop(i)
                ADList.pop(i)

    if TypeOf == 'tracer':
        return [IDs, Xs, Ys]
    elif TypeOf == 'individual':
        return [Type, Xs, Ys, IDs, ID, Vals, GrowthList, MaintList,
                N_RList, P_RList, C_RList, DispList, ADList]
    elif TypeOf == 'resource':
        return [Type, Xs, Ys, IDs, ID, Vals]




def maintenance(Sp_IDs, Xs, Ys, MD, MFD, RPD, EnvD, IDs, Qs, GrowthList,
        MaintList, N_RList, P_RList, C_RList, DispList, ADList):

    if Sp_IDs == []:
        return [Sp_IDs, Xs, Ys, IDs, Qs, GrowthList, MaintList, N_RList,
                P_RList, C_RList, DispList, ADList]

    n = len(IDs)
    for j in range(n):

        i = randint(0, len(IDs)-1)

        val = Qs[i]
        val[0] -= MaintList[i] # maintanence influenced by species id
        val[1] -= MaintList[i]
        val[2] -= MaintList[i]

        if min(val) <= MaintList[i]*0.00001:   # starved

            Qs.pop(i)
            Sp_IDs.pop(i)
            IDs.pop(i)
            Xs.pop(i)
            Ys.pop(i)
            GrowthList.pop(i)
            MaintList.pop(i)
            N_RList.pop(i)
            P_RList.pop(i)
            C_RList.pop(i)
            DispList.pop(i)
            ADList.pop(i)

        else: Qs[i] = val

    return [Sp_IDs, Xs, Ys, IDs, Qs, GrowthList, MaintList, N_RList,
            P_RList, C_RList, DispList, ADList]





def transition(Sp_IDs, IDs, Qs, GrowthList, MaintList, MFD, RPD, ADList):

    if Sp_IDs == []:
        return [Sp_IDs, IDs, Qs, GrowthList, MaintList, ADList]

    n = len(IDs)
    for j in range(n):

        i = randint(0, len(IDs)-1)
        spid = Sp_IDs[i]
        state = ADList[i]

        mfd = MFD[spid]

        if state == 'd':
            #continue
            x = np.random.binomial(1, RPD[spid]) # make this probability a randomly chosen variable
            if x == 1:

                ADList[i] = 'a'
                MaintList[i] = mfd*MaintList[i]

        if state == 'a':
            #continue
            val = Qs[i]
            if max(val) <= MaintList[i]*mfd:  # go dormant

                MaintList[i] = MaintList[i]/mfd # make this a randomly chosen variable
                ADList[i] = 'd'

    return [Sp_IDs, IDs, Qs, GrowthList, MaintList, ADList]





def consume(field, R_Types, R_Vals, R_IDs, R_ID, RXs, RYs, Sp_IDs, Qs, I_IDs, I_ID,
    IXs, IYs, w, h, GD, N_RD, P_RD, C_RD, DispD, GrowthList, MaintList, N_RList,
    P_RList, C_RList, DispList, ADList, enz_action = True):


    if Sp_IDs == []:
        return [Sp_IDs, I_IDs, Qs, GrowthList, MaintList, ADList]

    n = len(I_IDs)
    for i in range(n):

        x1 = IXs[i]
        y1 = IYs[i]

        r = len(R_IDs)
        Try = r
        ct = 0

        while ct < Try:
            ct += 1

            j = randint(0, r-1)
            x2 = RXs[j]
            y2 = RYs[j]

            dist = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

            if dist <= min(Qs[i]) * 1000:
                state = ADList[i]
                if state == 'd':
                    ADList[i] == 'a'

                # The food
                R_val = R_Vals[j]
                R_type = R_Types[j]

                rtype = list(R_type)
                R = rtype.pop(0)
                rnum = int(''.join(rtype))

                # The individual's cell quota
                Q = Qs[i]
                QN = Q[0]
                QP = Q[1]
                QC = Q[2]

                # the species
                sp = Sp_IDs[i]
                mu = GD[sp]

                Q = 0.0
                efficiency = 0.0

                if R == 'N':
                    efficiency = N_RList[i][rnum]
                    Q = QN

                if R == 'P':
                    efficiency = P_RList[i][rnum]
            	    Q = QP

                if R == 'C':
                    efficiency = C_RList[i][rnum]
                    Q = QC

                mu = mu * efficiency

                if R_val > (mu * Q): # Increase cell quota
                    R_val = R_val - (mu * Q)
                    Q += (mu * Q)

                else:
                    Q += R_val
                    R_val = 0.0

                if Q > 1.0:
                    R_val = Q - 1.0
                    Q = 1.0
                    R_Vals[j] = R_val


                if R_val <= 0.0:
                    R_Vals.pop(j)
                    R_Types.pop(j)
                    R_IDs.pop(j)
                    RXs.pop(j)
                    RYs.pop(j)

                if Q < 0.0:
                    print Q, QN, QP, QC
                    sys.exit()

                if R == 'N':
                    Qs[i] = [Q, QP, QC]
                if R == 'P':
                    Qs[i] = [QN, Q, QC]
                if R == 'C':
                    Qs[i] = [QN, QP, Q]

    return [R_Types, R_Vals, R_IDs, R_ID, RXs, RYs, Sp_IDs, Qs, I_IDs, I_ID, IXs,
        IYs, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList, ADList]






def consume_boxes(field, R_Types, R_Vals, R_IDs, R_ID, R_Xs, R_Ys, Sp_IDs,
        Qs, I_IDs, I_ID, I_Xs, I_Ys, w, h, GD, N_RD, P_RD, C_RD, DispD,
        GrowthList, MaintList, N_RList, P_RList, C_RList, DispList, ADList, enz_action = True):

    if not len(R_Types) or not len(Sp_IDs):
        List = [R_Types, R_Vals, R_IDs, R_ID, R_Xs]
        List += [R_Ys, Sp_IDs, Qs, I_IDs, I_ID]
        List += [I_Xs, I_Ys, GrowthList, MaintList, N_RList,
                P_RList, C_RList, DispList, ADList]
        return List

    I_Boxes = [list([]) for _ in xrange(w*h)]
    R_Boxes = [list([]) for _ in xrange(w*h)]

    index = 0
    for i, val in enumerate(I_IDs):
        rX = int(round(I_Xs[i]))
        rY = int(round(I_Ys[i]))

        index = int(round(rX + (rY * w)))

        if index > len(I_Boxes) - 1:
            index = len(I_Boxes) - 1
        elif index < 0:
            index = 0

        I_Boxes[index].append(val)

    index = 0
    for i, val in enumerate(R_IDs):

        rX = int(round(R_Xs[i]))
        rY = int(round(R_Ys[i]))
        index = int(round(rX + (rY * w)))

        if index > len(R_Boxes) - 1:
            index = len(R_Boxes) - 1
        elif index < 0:
            index = 0

        R_Boxes[index].append(val)


    for i, box in enumerate(I_Boxes):
        if not len(box): continue

        R_Box = R_Boxes[i]
        EnzVal = field[i] # enzymatic value

        for ind in box: # The individuals
            if not len(R_Box): break

            R_ID = choice(R_Box)
            boxI_ex = R_Box.index(R_ID)

            # The food
            j = R_IDs.index(R_ID)
            R_val = R_Vals[j]
            R_type = R_Types[j]

            rtype = list(R_type)
            R = rtype.pop(0)
            rnum = int(''.join(rtype))

            # The Individual
            ID = I_IDs.index(ind)

            state = ADList[ID]
            if state == 'd':
                ADList[ID] == 'a'

            # The individual's cell quota

            Q = Qs[ID]
            QN = Q[0]
            QP = Q[1]
            QC = Q[2]

            # the species
            sp = Sp_IDs[ID]
            mu = GD[sp]

            Q = 0.0
            efficiency = 0.0

            if R == 'N':
                efficiency = N_RList[ID][rnum]
                Q = QN

            if R == 'P':
                efficiency = P_RList[ID][rnum]
                Q = QP

            if R == 'C':
                efficiency = C_RList[ID][rnum]
                Q = QC

            mu = mu * efficiency

            if R_val > (mu * Q): # Increase cell quota
                R_val = R_val - (mu * Q)
                Q += (mu * Q)

            else:
                Q += R_val
                R_val = 0.0

            if Q > 1.0:
                R_val = Q - 1.0
                Q = 1.0
                R_Vals[j] = R_val


            if R_val <= 0.0:
                R_Box.pop(boxI_ex)
                R_Vals.pop(j)
                R_Types.pop(j)
                R_IDs.pop(j)
                R_Xs.pop(j)
                R_Ys.pop(j)


            if Q < 0.0:
                print Q, QN, QP, QC
                sys.exit()

            if R == 'N':
                Qs[ID] = [Q, QP, QC]
            if R == 'P':
                Qs[ID] = [QN, Q, QC]
            if R == 'C':
                Qs[ID] = [QN, QP, Q]


    return [R_Types, R_Vals, R_IDs, R_ID, R_Xs, R_Ys, Sp_IDs,
            Qs, I_IDs, I_ID, I_Xs, I_Ys, GrowthList, MaintList, N_RList,
            P_RList, C_RList, DispList, ADList]



def reproduce(spec, Sp_IDs, Qs, IDs, ID, Xs, Ys, w, h, GD, DispD,
        N_RD, P_RD, C_RD, MD, MFD, RPD, EnvD, envGs, nN, nP, nC, GList, MList,
        NList, PList, CList, DList, ADList):

    if Sp_IDs == []:
        return [Sp_IDs, Qs, IDs, ID, Xs, Ys, GD, DispD, GList, MList,
                NList, PList, CList, DList, ADList]

    repro = 'fission'
    if repro == 'fission':

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
                r1,r2,r3 = Qs[i]
                r1 -= MList[i]*r1
                r2 -= MList[i]*r2
                r3 -= MList[i]*r3
                Qs[i] = [r1, r2, r3]

                pg = []
                sp_opts = EnvD[spID]

                for g, opt in enumerate(sp_opts):

                    x, y = envGs[g]
                    pg.append(1 - (abs(X - x)/max([X,x])))
                    pg.append(1 - (abs(Y - y)/max([Y,y])))


                if np.mean(pg) > 1 or np.mean(pg) < 0:
                    print pg
                    sys.exit()

                p = np.mean(pg)
                p = np.random.binomial(1, p)
                if p == 1: # the environment is suitable for reproduction

                    QN = Q[0]
                    QP = Q[1]
                    QC = Q[2]

                    Qs[i] = [QN/2.0, QP/2.0, QC/2.0]
                    Qs.append([QN/2.0, QP/2.0, QC/2.0])

                    ID += 1
                    IDs.append(ID)

                    p = np.random.binomial(1, spec)
                    p = 0
                    if p == 1:

                        # speciate
                        #t = time.clock()
                        spID_new = str(spID)# +' '+ str(t)

                        # new species growth rate
                        p = np.random.binomial(1, 0.25)
                        if p == 1: GD[spID_new] = np.random.uniform(0.5, 1.0)
                        else: GD[spID_new] = GD[spID]

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

                        # new species resource use efficiencies
                        # Nitrogen
                        p = np.random.binomial(1, 0.25)
                        if p == 1: N_RD[spID_new] = np.random.uniform(0.01, 1.0, nN)
                        else: N_RD[spID_new] = N_RD[spID]

                        # Phosphorus
                        p = np.random.binomial(1, 0.25)
                        if p == 1: P_RD[spID_new] = np.random.uniform(0.01, 1.0, nP)
                        else: P_RD[spID_new] = P_RD[spID]

                        # Carbon
                        p = np.random.binomial(1, 0.25)
                        if p == 1: C_RD[spID_new] = np.random.uniform(0.01, 1.0, nC)
                        else: C_RD[spID_new] = C_RD[spID]

                        spID = spID_new

                    means = GD[spID]
                    i = GetIndParam(means)
                    GList.append(i)

                    means = MD[spID]
                    i = GetIndParam(means)
                    MList.append(i)

                    means = N_RD[spID]
                    i = GetIndParam(means)
                    NList.append(i)

                    means = P_RD[spID]
                    i = GetIndParam(means)
                    PList.append(i)

                    means = C_RD[spID]
                    i = GetIndParam(means)
                    CList.append(i)

                    means = DispD[spID]
                    i = GetIndParam(means)
                    DList.append(i)

                    Sp_IDs.append(spID)

                    ADList.append('a')

                    newX = float(np.random.uniform(X-0.1, X, 1))
                    if limit > newX: newX = 0
                    if newX > w - limit: newX = w - limit
                    Xs.append(newX)

                    newY = float(np.random.uniform(Y-0.1, Y+0.1, 1))
                    if limit > newY: newY = 0
                    elif newY > h: newY = h - limit
                    Ys.append(newY)


    listlen = [len(Sp_IDs), len(Qs), len(IDs), len(Xs), len(Ys), len(GList), len(MList), len(DList), len(ADList)]
    if min(listlen) != max(listlen):
        print 'In reproduce'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()

    return [Sp_IDs, Qs, IDs, ID, Xs, Ys, GD, DispD, GList, MList,
                NList, PList, CList, DList, ADList]




def chemotaxis(repro, spec, Sp_IDs, Qs, IDs, ID, Xs, Ys,  w, h, GD, DispD,
        N_RD, P_RD, C_RD, MD, MFD, RPD, EnvD, envGs, nN, nP, nC, GList, MList,
        NList, PList, CList, DList, ADList):

    if Sp_IDs == []:
        return [Sp_IDs, Qs, IDs, ID, Xs, Ys, GD, DispD, GList,
                    MList, NList, PList, CList, DList, ADList]

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

                r1,r2,r3 = Qs[i]
                r1 -= MList[i]*dist*r1
                r2 -= MList[i]*dist*r2
                r3 -= MList[i]*dist*r3
                Qs[i] = [r1, r2, r3]

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


    return [Sp_IDs, Qs, IDs, ID, Xs, Ys, GD, DispD, GList, MList,
                NList, PList, CList, DList, ADList]



def density_forage(RVals, RX, RY, repro, spec, Sp_IDs, Qs, IDs, ID, Xs, Ys,  w, h, GD, DispD,
        N_RD, P_RD, C_RD, MD, MFD, RPD, EnvD, envGs, nN, nP, nC, GList, MList,
        NList, PList, CList, DList, ADList):

    if Sp_IDs == []:
        return [Sp_IDs, Qs, IDs, ID, Xs, Ys, GD, DispD, GList,
                    MList, NList, PList, CList, DList, ADList]

    # Locate resource density
    # 1.) Get mean X and Y values for resources
    avgX = np.mean(RX)
    avgY = np.mean(RY)

    n = len(IDs)
    for j in range(n):

        i = randint(0, len(IDs)-1)

        state = ADList[i]
        if state == 'd':
            continue

        spID = Sp_IDs[i]
        X = Xs[i]
        Y = Ys[i]

        dist = DispD[spID]

        # A cost for active dispersal
        r1,r2,r3 = Qs[i]
        r1 -= MList[i]*dist*r1
        r2 -= MList[i]*dist*r2
        r3 -= MList[i]*dist*r3
        Qs[i] = [r1, r2, r3]


        if X > avgX:
            X -= dist

        elif X < avgX:
            X += dist

        if Y > avgY:
            Y -= dist

        elif Y < avgY:
            Y += dist

        if X > w: X = w
        elif X < 0: X = 0
        if Y > h: Y = h
        elif Y < 0: Y = 0

        Xs[i] = X
        Ys[i] = Y

    return [Sp_IDs, Qs, IDs, ID, Xs, Ys, GD, DispD, GList, MList,
                NList, PList, CList, DList, ADList]


def nearest_forage(RVals, RX, RY, repro, spec, Sp_IDs, Qs, IDs, ID, Xs, Ys,  w, h, GD, DispD,
        N_RD, P_RD, C_RD, MD, MFD, RPD, EnvD, envGs, nN, nP, nC, GList, MList,
        NList, PList, CList, DList, ADList):

    if Sp_IDs == []:
        return [Sp_IDs, Qs, IDs, ID, Xs, Ys, GD, DispD, GList,
                    MList, NList, PList, CList, DList, ADList]

    n = len(IDs)
    n = min([100, n])
    r = len(RVals)

    for j in range(n):
        i = randint(0, len(IDs)-1)

        state = ADList[i]
        if state == 'd':
            continue

        x1 = Xs[i]
        y1 = Ys[i]

        MinDist = 10000

        rx = 0
        ry = 0

        for j in range(r):

            x2 = RX[j]
            y2 = RY[j]

            dist = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)
            if dist < MinDist:
                MinDist = dist
                rx = x2
                ry = y2

        spID = Sp_IDs[i]
        dist = DispD[spID]

        # A cost for active dispersal
        r1,r2,r3 = Qs[i]
        r1 -= MList[i]*dist*r1
        r2 -= MList[i]*dist*r2
        r3 -= MList[i]*dist*r3
        Qs[i] = [r1, r2, r3]


        if x1 > rx:
            x1 -= dist

        elif x1 < rx:
            x1 += dist

        if y1 > ry:
            y1 -= dist

        elif y1 < ry:
            y1 += dist

        if x1 > w: x1 = w
        elif x1 < 0: x1 = 0
        if y1 > h: y1 = h
        elif y1 < 0: y1 = 0

        Xs[i] = x1
        Ys[i] = y1

    return [Sp_IDs, Qs, IDs, ID, Xs, Ys, GD, DispD, GList, MList,
                NList, PList, CList, DList, ADList]
