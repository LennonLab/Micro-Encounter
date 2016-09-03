from __future__ import division
import statsmodels.tsa.stattools as sta
import linecache
from random import choice
from math import isnan
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
ComplexityLevels = [choice([1,2,3]), choice([1,2,3,4]), choice([1,2,3])]
Mu, Maint, ct, IndID, RID, N, T, R, PRODI, PRODQ, numD = [0]*11
ADList, ADs, AVG_DIST, SpecDisp, SpecMaint, SpecGrowth, Deadlist = [list([]) for _ in xrange(7)]
IndLists, ResLists, Gs, Ms, Qs, Ds, Rs, PRODIs, Ns, RDENs, MUs, MAINTs, encList, Ragg, Iagg = [list([]) for _ in xrange(15)]

# GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict
IndDicts = [{}, {}, {}, {}, {}, {}]
ResLists = [[], [], [], []]
IndLists = [[], [], [], [], [], []]

#######################  Other variables  ######################################
num_sims, LowerLimit, sim, p, ct = 100000, 30, 1, 0, 0
BurnIn = 'not done'
Nlim = 1000

#######################  Main Simulation Loop  #################################

while sim < num_sims:

    numDead, encounters = 0, 0
    ct += 1
    
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
    numD = ADList.count('d')


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

        #print sim, ct, N
        PRODIs.append(PRODI)
        Rvals, RX, RY, RIDs = ResLists

        RDENs.append(len(Rvals)/(height*width))
        Rs.append(len(RX))
        encList.append(encounters)

        if N > 1:
            SpeciesIDs, IndX, IndY, IndIDs, CellQuotas, ADList = IndLists
            Iagg.append(spatial.morisitas(IndX, IndY, width, height))

            if R > 1:
                q = min([20, R])
                AVG_DIST.append(spatial.nearest_neighbor(IndX, RX, IndY, RY, q))
                Ragg.append(spatial.morisitas(RX, RY, width, height))

            GrowthDict, MaintDict, MainFactorDict, RPFDict, DispDict, TrophicDict = IndDicts

            SpecDisp.append(mean(DispDict.values()))
            SpecMaint.append(mean(MaintDict.values()))
            SpecGrowth.append(mean(GrowthDict.values()))

            Gs.append(bide.per_capita(GrowthDict, SpeciesIDs))
            Ms.append(bide.per_capita(MaintDict, SpeciesIDs))
            Ds.append(bide.per_capita(DispDict, SpeciesIDs))
            Qs.append(mean(CellQuotas))

            numD = ADList.count('d')
            ADs.append(numD/len(ADList))
            Deadlist.append(numDead)

        if len(Ns) >= 20 or N > Nlim:

            if len(Ns) >= 20:

                print '%4s' % sim, ' r:','%4s' %  r, ' R:','%4s' % int(round(mean(Rs))), ' N:','%5s' % int(round(mean(Ns))), \
                ' Dormant:', '%5s' % round(mean(ADs),3), ' Encounters:','%5s' % round(mean(encList),2), '   Spatial:', ComplexityLevels[0],\
                ' Trophic:', ComplexityLevels[1], ' Resource:', ComplexityLevels[2]

                outlist = [sim, mean(PRODIs), var(PRODIs), r, gmax, maintmax, dmax, seedCom, \
                width-0.2, height, mean(Ns), var(Ns), m, mean(RDENs), var(RDENs), mean(Rs), var(Rs), \
                mean(Gs), var(Gs), mean(Ms), var(Ms), mean(Ds), var(Ds), \
                mean(SpecGrowth), var(SpecGrowth), mean(SpecDisp), var(SpecDisp), \
                mean(SpecMaint),  var(SpecMaint), mean(AVG_DIST), var(AVG_DIST), \
                mean(ADs), var(ADs), ComplexityLevels[1], ComplexityLevels[2], \
                ComplexityLevels[0], mean(encList), var(encList), std, mean(Iagg), \
                var(Iagg), mean(Ragg), var(Ragg), mean(Deadlist), var(Deadlist), ct]

                outlist = str(outlist).strip('[]')
                outlist = str(outlist).strip('')

                OUT = open(mydir + '/results/simulated_data/SimData.csv', 'a')
                print>>OUT, outlist
                OUT.close()

            ComplexityLevels = [choice([1,2,3]), choice([1,2,3,4]), choice([1,2,3])]
            Mu, Maint, ct, IndID, RID, N, T, R, PRODI, PRODQ, numD = [0]*11
            ADList, ADs, AVG_DIST, SpecDisp, SpecMaint, SpecGrowth, Deadlist = [list([]) for _ in xrange(7)]
            IndLists, ResLists, Gs, Ms, Qs, Ds, Rs, PRODIs, Ns, RDENs, MUs, MAINTs, encList, Ragg, Iagg = [list([]) for _ in xrange(15)]
            width, height, seedCom, m, r, gmax, maintmax, dmax, pmax, mmax, std = rp.get_rand_params()

            IndDicts = [{}, {}, {}, {}, {}, {}]
            ResLists = [[], [], [], []]
            IndLists = [[], [], [], [], [], []]

            p = 0
            BurnIn = 'not done'
            sim += 1
            ct = 0
