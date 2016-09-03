from __future__ import division
import statsmodels.tsa.stattools as sta
import linecache
from random import choice
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
sys.path.append(mydir + "/model/col_labels")

labels = linecache.getline(mydir + '/model/col_labels/labels.txt', 1)
with open(mydir + '/results/simulated_data/SimData.csv', 'w+') as text_file:
    text_file.write(labels)


#########################  Randomly chosen variables  ##########################

# Complexity levels: spatial, trophic, resource
ComplexityLevels = [choice([1,2,3]), choice([1,2,3,4]), choice([1,2,3])]

#Randomly chosen values for: width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std
params = rp.get_rand_params()
width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std = params

#######################  Lists & Dictionaries  #################################
ComplexityLevels = [choice(['whitenoise', 'aggregated_resources', 'chemo', 'fluid']),
                    choice(['null_comm', 'crossfeeding', 'scavenging', 'pred-prey', 'mutualism']),
                    choice(['monoculture', 'polyculture', 'lockandkey'])]

ct, IndID, RID, N = [0]*4
ADList, ADs, AVG_DIST, SpecDisp, SpecMaint, SpecGrowth, Deadlist = [list([]) for _ in xrange(7)]
IndLists, ResLists, Ns = [list([]) for _ in xrange(3)]

IndDicts = [{}, {}, {}, {}, {}, {}]
IndLists = [[], [], [], [], [], []]
ResLists = [[], [], [], []]
RDict = {}

if 'lockandkey' in ComplexityLevels[2] :
    RDict['dead'] = np.random.uniform(0.5, 1.0)
else:
    RDict['dead'] = np.random.uniform(1.0, 1.0)


#######################  Other variables  ######################################
num_sims, LowerLimit, sim, p, ct = 100000, 30, 1, 0, 0
BurnIn = 'not done'
Nlim = 1000

#######################  Main Simulation Loop  #################################

while sim < num_sims:

    numDead, encounters = 0, 0

    # Inflow of resources
    ResLists, RID = bide.ResIn(ResLists, RID, params, ct, ComplexityLevels)

    # Immigration
    IndLists, IndDicts, IndID = bide.immigration(IndLists, IndDicts, IndID, params, ct, ComplexityLevels)

    # Individual Dispersal
    if ComplexityLevels[0] < 3 and len(IndLists[0]) > 0: # spatial complexity
        IndList, IndDicts, IndID, ResList, RID, numDead = bide.dispersal(IndLists, IndDicts, IndID, ResLists, RID, params, ComplexityLevels, numDead)

    elif ComplexityLevels[0] == 3 and len(IndLists[0]) > 0: # spatial complexity
        ResList, RID, IndList, IndDicts, IndID, numDead = bide.chemotaxis(IndLists, IndDicts, IndID, ResLists, RID, params, ComplexityLevels, numDead)

    # Resource Dispersal
    if ComplexityLevels[0] == 1 and len(ResLists[0]) > 0: # spatial complexity
        ResLists, RID = bide.res_dispersal(ResLists, RID, params, ct, ComplexityLevels)

    p1 = len(IndList[0])

    # Consumption
    if len(IndLists[0]) > 0:
        ResLists, RID, IndLists, encounters = bide.consume(IndLists, IndDicts, IndID, ResLists, RID, params, ComplexityLevels)

    # Reproduction
    if len(IndLists[0]) > 0:
        IndLists, IndDicts, IndID, ResLists, RID, numDead = bide.reproduce(IndLists, IndDicts, IndID, ResLists, RID, params, ComplexityLevels, numDead)

    # Maintenance
    if len(IndLists[0]) > 0:
        IndLists, ResLists, RID, numDead = bide.maintenance(IndLists, IndDicts, ResLists, RID, ComplexityLevels, numDead)

    # Transition to or from dormancy
    if len(IndLists[0]) > 0:
        IndLists, ResLists, RID, numDead = bide.transition(IndLists, IndDicts, ResLists, RID, ComplexityLevels, numDead)

    p2 = len(IndLists[0])
    PRODI = p2 - p1

    N = len(IndLists[0])
    Ns.append(N)

    if N > Nlim:
        BurnIn = 'done'
        Ns = [Ns[-1]] # only keep the most recent N value

    if len(Ns) >= 50:

        if BurnIn == 'not done':
            AugmentedDickeyFuller = sta.adfuller(Ns)
            val, p = AugmentedDickeyFuller[0:2]

            if p >= 0.05:
                Ns.pop(0)

            elif p < 0.05 or isnan(p) == True:
                BurnIn = 'done'
                Ns = [Ns[-1]] # only keep the most recent N value

    if BurnIn == 'done':
        ct += 1

        Rvals, RX, RY, RIDs = ResLists
        SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList = IndLists
        GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = IndDicts
        Iagg = spatial.morisitas(IndX, IndY, width, height)

        if len(RX) > 1: Ragg = spatial.morisitas(RX, RY, width, height)

        outlist = [sim, ct, width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std, \
            ComplexityLevels[1], ComplexityLevels[2], ComplexityLevels[0],\
            mean(DispDict.values()), mean(MaintDict.values()), mean(GrowthDict.values()), \
            bide.per_capita(GrowthDict, SpeciesIDs), bide.per_capita(MaintDict, SpeciesIDs),\
            bide.per_capita(DispDict, SpeciesIDs), mean(CellQuotas), ADList.count('d'), numDead,\
            PRODI, N, len(Rvals)/(height*width), len(RX), encounters, Iagg, Ragg]

        outlist = str(outlist).strip('[]')
        outlist = str(outlist).strip('')
        OUT = open(mydir + '/results/simulated_data/SimData.csv', 'a')
        print>>OUT, outlist
        OUT.close()

        if ct > 100 or N > Nlim:
            if ct > 100:

                print '%4s' % sim, ' r:','%4s' %  r, ' R:','%4s' % len(RX), ' N:','%5s' % int(round(mean(Ns))), \
                ' Dormant:', '%5s' % round(ADList.countd('d')/N,3), ' Encounters:','%5s' % encounters,\
                '   Spatial:', ComplexityLevels[0], ' Trophic:', ComplexityLevels[1], ' Resource:', ComplexityLevels[2]

            ComplexityLevels = [choice(['whitenoise', 'aggregated_resources', 'chemo', 'fluid']),
                    choice(['null_comm', 'crossfeeding', 'scavenging', 'pred-prey', 'mutualism']),
                    choice(['monoculture', 'polyculture', 'lockandkey'])]

            ct, IndID, RID, N = [0]*5
            width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std = rp.get_rand_params()

            IndDicts = [{}, {}, {}, {}, {}, {}]
            IndLists = [[], [], [], [], [], []]
            ResLists = [[], [], [], []]
            ResDict = {}

            BurnIn = 'not done'
            sim += 1
