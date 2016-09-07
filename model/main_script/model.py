from __future__ import division
import statsmodels.tsa.stattools as sta
import linecache
from random import choice, randint, shuffle, sample
from math import isnan
import numpy as np
from numpy import mean, var
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

sys.path.append(mydir + "/model/col_labels")
labels = linecache.getline(mydir + '/model/col_labels/labels.txt', 1)
with open(mydir + '/results/simulated_data/SimData.csv', 'w+') as text_file:
    text_file.write(labels)

#########################  Randomly chosen variables  ##########################
ComplexityLevels = metrics.get_complexity_levels()
SC, TC, RC = ComplexityLevels

params = rp.get_rand_params()
width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std = params

#######################  Variables, Lists & Dictionaries  ######################
Ns = []
ResDict = {}
ct, IndID, RID, N = [0]*4
x = range(0, 8)
IndDicts = [{}, {}, {}, {}, {}, {}]
IndLists = [[], [], [], [], [], []]
ResLists = [[], [], [], [], []]

if 'lockandkey' in RC: 
    ResDict['dead'] = round(np.random.uniform(0.05, 0.1), 3)
elif 'simple' in RC: 
    ResDict['dead'] = 1.0

Nlim, sim, p, ct, RowID = 10000, 1, 0, 0, 0
BurnIn = 'not done'

#######################  Main Simulation Loop  #################################

while sim:
    numDead, encounters = 0, 0
    ct += 1
    RowID += 1

    shuffle(x)
    for xi in x:
        # Inflow of resources
        if xi == 0: ResLists, ResDict, RID = bide.ResIn(ResLists, ResDict, RID, params, ct, ComplexityLevels)

        # Immigration
        if xi == 1: IndLists, IndDicts, IndID = bide.immigration(IndLists, IndDicts, IndID, params, ct, ComplexityLevels)

        # Individual Dispersal
        if xi == 2 and 'nodispersal' not in SC: IndLists, IndDicts, IndID, ResList, ResDict, RID, numDead = bide.dispersal(IndLists, IndDicts, IndID, ResLists, ResDict, RID, params, ComplexityLevels, numDead)

        # Resource Dispersal
        if xi == 3 and 'static' not in SC: ResLists, ResDict, RID = bide.res_dispersal(ResLists, ResDict, RID, params, ct, ComplexityLevels)
    
        # Consumption
        if xi == 4: ResLists, ResDict, RID, IndLists, encounters, numDead = bide.consume(IndLists, IndDicts, IndID, ResLists, ResDict, RID, params, ComplexityLevels, numDead)

        # Reproduction
        if xi == 5: PRODI, IndLists, IndDicts, IndID, ResLists, ResDict, RID, numDead = bide.reproduce(IndLists, IndDicts, IndID, ResLists, ResDict, RID, params, ComplexityLevels, numDead)

        # Maintenance
        if xi == 6: IndLists, ResLists, ResDict, RID, numDead = bide.maintenance(IndLists, IndDicts, ResLists, ResDict, RID, ComplexityLevels, numDead)
        
        # Transition to or from dormancy
        if xi == 7: IndLists, ResLists, ResDict, RID, numDead = bide.transition(IndLists, IndDicts, ResLists, ResDict, RID, ComplexityLevels, numDead)

    N = len(IndLists[0])
    Ns.append(N)
    
    if len(Ns) >= 100 and BurnIn == 'not done':
        
            AugmentedDickeyFuller = sta.adfuller(Ns)
            val, p = AugmentedDickeyFuller[0:2]

            if p >= 0.05:
                Ns.pop(0)

            elif p < 0.05 or isnan(p) == True:
                BurnIn = 'done'
                Ns = [Ns[-1]] # only keep the most recent N value
                ct = 0
                
    if BurnIn == 'done':
        Rvals, Rtypes, RX, RY, RIDs = ResLists
        SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList = IndLists
        GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = IndDicts
        
        Iagg = spatial.morisitas(IndX, IndY, width, height)
        Ragg = spatial.morisitas(RX, RY, width, height)

        outlist = [RowID, sim, ct, width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std, \
            ComplexityLevels[0], ComplexityLevels[1], ComplexityLevels[2],\
            mean(DispDict.values()), mean(MaintDict.values()), mean(GrowthDict.values()), \
            metrics.per_capita(GrowthDict, SpeciesIDs), metrics.per_capita(MaintDict, SpeciesIDs),\
            metrics.per_capita(DispDict, SpeciesIDs), mean(CellQuotas), ADList.count('d'), numDead,\
            PRODI, N, len(RX), encounters, Iagg, Ragg]

        outlist = str(outlist).strip('[]')
        outlist = str(outlist).strip('')
        OUT = open(mydir + '/results/simulated_data/SimData.csv', 'a')
        print>>OUT, outlist
        OUT.close()

        if len(Ns) > 50 or N > Nlim:
            if len(Ns) > 50:
                print 'sim:',sim, ' N:',N, ' R:',len(ResLists[0]), ' %Dormant:',100*round(IndLists[-1].count('d')/N, 3), ' Encounters:',encounters, ' Prod:',PRODI, ' dead:',numDead
                    
            ComplexityLevels = metrics.get_complexity_levels()
            SC, TC, RC = ComplexityLevels
            
            ct, IndID, RID, N = [0]*4
            width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std = rp.get_rand_params()
    
            IndDicts = [{}, {}, {}, {}, {}, {}]
            IndLists = [[], [], [], [], [], []]
            ResLists = [[], [], [], [], []]
            if 'lockandkey' in RC: 
                ResDict['dead'] = round(np.random.uniform(0.1, 0.5), 3)
            elif 'simple' in RC: 
                ResDict['dead'] = 1.0

            BurnIn = 'not done'
            sim += 1