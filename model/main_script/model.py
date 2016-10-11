from __future__ import division
import statsmodels.tsa.stattools as sta
import linecache
from random import shuffle
from math import isnan
import numpy as np
from numpy import mean
import time
import sys
import os

mydir = os.path.expanduser("~/GitHub/Micro-Encounter")
sys.path.append(mydir + "/model/bide")
import bide
sys.path.append(mydir + "/model/randparams")
import randparams as rp
sys.path.append(mydir + "/model/spatial")
import spatial
sys.path.append(mydir + "/model/metrics")
import metrics

'''
sys.path.append(mydir + "/model/col_labels")
labels = linecache.getline(mydir + '/model/col_labels/labels.txt', 1)
with open(mydir + '/results/simulated_data/SimData.csv', 'w+') as text_file:
    text_file.write(labels)
'''

#########################  Randomly chosen variables  ##########################
ComplexityLevels = metrics.get_complexity_levels()
SC, TC, RC = ComplexityLevels

extremes = True
params = rp.get_rand_params(extremes)
width, height, length, seed, m, r, gmax, mmax, dmax, pmax, mfact, std = params

#######################  Variables, Lists & Dictionaries  ######################
Ns = []
IndDict, ResDict = {}, {}
ct, N = 0, 0
x = range(0, 8)
SpDicts = [{}, {}, {}, {}, {}, {}]
ResLists = [[], [], [], [], [], []]

if 'lockandkey' in RC:
    ResDict['dead'] = np.random.uniform(0.1, 0.5)
elif 'simple' in RC:
    ResDict['dead'] = 1.0

numSims, Nlim, sim, p, ct, RowID, PRODI = 10**6, 1000, 0, 0, 0, 0, 0
BurnIn = 'not done'

#######################  Main Simulation Loop  #################################
t0 = time.clock()
while sim < numSims:
    numDead, encounters, res_in = 0, 0, 0
    ct += 1
    RowID += 1

    shuffle(x)
    t2 = float()
    t1 = time.clock()
    for xi in x:

        # Inflow of resources
        if xi == 0: ResLists, ResDict, res_in = bide.ResIn(ResLists, ResDict, params, ct, ComplexityLevels)

        # Immigration
        elif xi == 1: IndDict, SpDicts = bide.immigration(IndDict, SpDicts, params, ct, ComplexityLevels)

        # Individual Dispersal
        elif xi == 2 and '-none-' not in SC: IndDict, ResList, ResDict, numDead = bide.dispersal(IndDict, SpDicts, ResLists, ResDict, params, ComplexityLevels, numDead)

        # Resource Dispersal
        elif xi == 3: ResLists = bide.res_dispersal(ResLists, params, ct, ComplexityLevels)

        # Consumption
        elif xi == 4: ResLists, ResDict, IndDict, encounters, numDead = bide.consume(IndDict, SpDicts, ResLists, ResDict, params, ComplexityLevels, numDead)

        # Reproduction
        elif xi == 5: PRODI, IndDicts, ResLists, ResDict, numDead = bide.reproduce(IndDict, SpDicts, ResLists, ResDict, params, ComplexityLevels, numDead)

        # Maintenance
        elif xi == 6: IndDict, ResLists, ResDict, numDead = bide.maintenance(IndDict, SpDicts, ResLists, ResDict, ComplexityLevels, numDead)

        # Transition to or from dormancy
        elif xi == 7: IndDict, ResLists, ResDict, numDead = bide.transition(IndDict, SpDicts, ResLists, ResDict, ComplexityLevels, numDead)

        t2 = time.clock() - t1
        if t2 >= 0.12: break

    #if t2 < 0.12: time.sleep(0.12 - t2)
    N = len(IndDict.keys())
    Ns.append(N)
    Rvs = ResLists[0]
    Rvs.sort()
    R = len(Rvs)

    if ct%100 == 0:
        print 'sim:', sim, 'ct:', ct, ' N:', N, ' R:', R, 'prodi:', PRODI, ' dead:', numDead, 'encounters:', encounters,' ', 'pmax:', round(pmax,2), 'mfact:', mfact

    BurnIn = 'done'

    if len(Ns) >= 500 and BurnIn == 'not done':
            AugmentedDickeyFuller = sta.adfuller(Ns)
            val, p = AugmentedDickeyFuller[0:2]

            if p >= 0.05: Ns.pop(0)

            elif p < 0.05 or isnan(p) == True:
                BurnIn = 'done'
                Ns = [Ns[-1]] # only keep the most recent N value
                ct = 0

    if BurnIn == 'done' and ct%10 == 0:
        Rvals, Rtypes, RX, RY, RZ, RIDs = ResLists

        IndIDs = IndDict.keys()
        SpeciesIDs, IndX, IndY, IndZ, ADList, CellQuotas = [], [], [], [], [], []

        [SpeciesIDs.append(IndDict[i]['species']) for i in IndIDs]
        [IndX.append(IndDict[i]['x']) for i in IndIDs]
        [IndY.append(IndDict[i]['y']) for i in IndIDs]
        [IndZ.append(IndDict[i]['z']) for i in IndIDs]
        [CellQuotas.append(IndDict[i]['quota']) for i in IndIDs]
        [ADList.append(IndDict[i]['state']) for i in IndIDs]

        GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = SpDicts
        Iagg = 1 #spatial.morisitas(IndX, IndY, width, height, length)
        Ragg = 1 #spatial.morisitas(RX, RY, width, height, length)

        outlist = [RowID, sim, ct, width, height, length, seed, m, r, gmax, mmax, dmax, pmax, mfact, std, \
            res_in, ComplexityLevels[0], ComplexityLevels[1], ComplexityLevels[2],\
            mean(DispDict.values()), mean(MaintDict.values()), mean(GrowthDict.values()), \
            metrics.per_capita(GrowthDict, SpeciesIDs), metrics.per_capita(MaintDict, SpeciesIDs),\
            metrics.per_capita(DispDict, SpeciesIDs), mean(CellQuotas), ADList.count('dormant'), numDead,\
            PRODI, N, len(RX), encounters, Iagg, Ragg]

        outlist = str(outlist).strip('[]')
        outlist = str(outlist).strip('')
        OUT = open(mydir + '/results/simulated_data/SimData.csv', 'a')
        print>>OUT, outlist
        OUT.close()

        limlist = [N, R]
        if len(Ns) > 2000 or min(limlist) > Nlim:

            if len(Ns) > 2000:
                N = int(round(np.mean(Ns)))
                if N == 0:
                    N = 1

                print 'sim:',sim, ' N:',N, ' R:',len(ResLists[0]), ' %Dormant:',100*round(ADList.count('dormant')/N, 3), ' Encounters:',encounters, ' Prod:',PRODI, ' dead:',numDead, ' SC:', SC

            ComplexityLevels = metrics.get_complexity_levels()
            SC, TC, RC = ComplexityLevels

            ct, N = 0, 0
            extremes = True
            params = rp.get_rand_params(extremes)
            width, height, length, seed, m, r, gmax, mmax, dmax, pmax, mfact, std = params
            #print params

            SpDicts = [{}, {}, {}, {}, {}, {}]
            IndDict, ResDict = {}, {}

            ResLists = [[], [], [], [], [], []]
            Ns = []

            if '-lockandkey-' in RC:
                ResDict['dead'] = np.random.uniform(0.1, 0.5)
            elif '-simple-' in RC:
                ResDict['dead'] = 1.0

            BurnIn = 'not done'
            sim += 1
