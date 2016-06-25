from __future__ import division
import matplotlib.animation as animation
import matplotlib.pyplot as plt

import statsmodels.tsa.stattools as sta
from random import choice, randint
from math import isnan

import numpy as np
from numpy import mean, var
import sys
import os

mydir = os.path.expanduser("~/")
sys.path.append(mydir + "GitHub/Micro-Encounter/model/bide")
import bide
sys.path.append(mydir + "GitHub/Micro-Encounter/model/randparams")
import randparams as rp
sys.path.append(mydir + "GitHub/Micro-Encounter/model/spatial")
import spatial

GenPath = mydir + 'GitHub/Micro-Encounter/results/simulated_data/'

OUT1 = open(GenPath + 'SimData.csv','w')
print>>OUT1,'RowID,MeanIndProduction,VarIndProduction,ResInflow,MaxGrowthRate,MaxMetMaint,MaxActiveDispersal,LogseriesA,StartingSeed,\
width,height,MeanTotalAbundance,VarTotalAbundance,ImmigrationRate,MeanResourceConcentration,VarResourceConcentration,\
MeanResourceParticles,VarResourceParticles,Speciation,MeanPerCapitaGrowth,VarPerCapitaGrowth,MeanPerCapitaMaint,\
VarPerCapitaMaint,MeanPerCapitaActiveDispersal,VarPerCapitaActiveDispersal,MeanSpecGrowth,VarSpecGrowth,\
MeanSpecDisp,VarSpecDisp,MeanSpecMaint,VarSpecMaint,MeanAvgDist,VarAvgDist,MeanDormFreq,VarDormFreq,\
TrophicComplexityLevel,ResourceComplexityLevel,SpatialComplexityLevel,MeanEncounter,VarEncounter,\
IncomingResAgg,MeanIndAgg,VarIndAgg,MeanResAgg,VarResAgg,MeanDeaths,VarDeaths,RunTime'
OUT1.close()

def nextFrame(arg):

    """ Function called for each successive animation frame; arg is the frame number """

    global mmax, pmax, ADs, ADList, AVG_DIST, SpecDisp, SpecMaint, SpecGrowth
    global fixed, p, BurnIn, t, num_sims, width, height, Rates, GrowthDict, RD
    global DispDict, MaintDict, gmax, dmax, maintmax, IndIDs, Qs, EVList
    global IndID, IndX, IndY, Ind_scatImage, SpeciesIDs, TLList
    global RX, RY, RID, RIDs, RVals, EnvD, resource_scatImage, Mu, Maint
    global speciation, seedCom, m, r, sim
    global N, ct, RDens, RDiv, RRich, T, R, LowerLimit, prod_i, prod_q, alpha
    global Ts, Rs, PRODIs, Ns, RDENs, RDIVs, RRICHs, GrowthList, MaintList, RList
    global MUs, MAINTs, PRODNs, PRODPs, PRODCs, Gs, Ms, Ds
    global DispList, envgrads, MainFactorDict, RPFDict, enzyme_field, u0

    global TrophicComplexityLevel, SpatialComplexityLevel, encList, std
    global ResourceComplexityLevel, BiologicalComplexityLevel
    global Ragg, Iagg, static, maxgen, Deadlist

    numDead = 0
    ct += 1
    #print ct
    plot_system = 'no'

    listlen = [len(SpeciesIDs), len(Qs), len(IndIDs), len(IndX), len(IndY), len(GrowthList), len(MaintList), len(DispList), len(ADList), len(EVList), len(TLList)]
    if min(listlen) != max(listlen):
        print ct, 'In model.py'
        print 'min(listlen) != max(listlen)'
        print listlen
        sys.exit()

    encounters = 0

    # Inflow of resources
    RList, RVals, RX, RY,  RIDs, RID = bide.ResIn(std, ct, RList, RVals, RX, RY,  RID,\
    RIDs, r, width, height, u0, TrophicComplexityLevel, SpatialComplexityLevel, \
    ResourceComplexityLevel, BiologicalComplexityLevel)

    # Immigration
    SpeciesIDs, IndX, IndY,  MaintDict, MainFactorDict, RPFDict, EnvD, GrowthDict,\
    DispDict, IndIDs, IndID, Qs, RD, GrowthList, MaintList, DispList, ADList, EVList,\
    TLList = bide.immigration(maxgen, std, mmax, pmax, dmax, gmax, maintmax, seedCom, m, \
    SpeciesIDs, IndX, IndY, width, height, MaintDict, MainFactorDict, RPFDict, EnvD, envgrads,\
    GrowthDict, DispDict, IndIDs, IndID, Qs, RD, u0, alpha, GrowthList, MaintList, \
    DispList, ADList, EVList, TLList, ct, TrophicComplexityLevel, \
    SpatialComplexityLevel, ResourceComplexityLevel, BiologicalComplexityLevel)

    # Dispersal
    if SpatialComplexityLevel < 3:
        SpeciesIDs, Qs, IndIDs, ID, IndX, IndY, GrowthDict, DispDict, GrowthList, \
        MaintList, DispList, ADList, EVList, TLList, RList, RVals, RX, RY, RIDs, RID, numDead = bide.dispersal(speciation, SpeciesIDs,\
	Qs, IndIDs, IndID, IndX, IndY, width, height, GrowthDict, \
        DispDict, RD, MaintDict, MainFactorDict, RPFDict, EnvD, envgrads, \
        GrowthList, MaintList, DispList, ADList, EVList, TLList, RList, RVals, RX, RY, RIDs, RID, TrophicComplexityLevel, \
        SpatialComplexityLevel, ResourceComplexityLevel, BiologicalComplexityLevel, numDead)

    elif SpatialComplexityLevel == 3:
        RList, RVals, RIDs, RID, RX, RY, SpeciesIDs, Qs, IndIDs, IndID, IndX, IndY, width, height, GD, RD, DispD, GrowthList,\
        MaintList, DispList, ADList, TLList, EVList, numDead = bide.chemotaxis(RList, RVals, RIDs, RID, RX, RY, SpeciesIDs, Qs, IndIDs, IndID,\
        IndX, IndY, width, height, GrowthDict, RD, DispDict, GrowthList, MaintList, DispList, ADList, TLList, EVList,\
        TrophicComplexityLevel, SpatialComplexityLevel, ResourceComplexityLevel, BiologicalComplexityLevel, numDead)

    # Resource Dispersal
    if SpatialComplexityLevel == 1:
        RList, RVals, RX, RY, RID, RIDs = bide.res_dispersal(ct, RList, RVals, RX, RY, RID, RIDs, r,\
        width, height, u0, TrophicComplexityLevel, SpatialComplexityLevel, \
        ResourceComplexityLevel, BiologicalComplexityLevel)

    PRODI = 0
    p1 = len(Qs)

    # Modify enzyme field
    # enzyme_field = field.EnzymeField(enzyme_field, IndX, IndY, ADList, Qs, width)

    # Consume
    RList, RVals, RIDs, RID, RX, RY, SpeciesIDs, Qs, encounters = bide.consume(enzyme_field, \
    RList, RVals, RIDs, RID, RX, RY, SpeciesIDs, Qs, IndIDs, IndID, IndX, IndY, \
    width, height, GrowthDict, RD, DispDict, GrowthList, MaintList, DispList, ADList, \
    TLList, EVList, TrophicComplexityLevel, SpatialComplexityLevel, \
    ResourceComplexityLevel, BiologicalComplexityLevel)

    # Reproduction
    SpeciesIDs, Qs, IndIDs, ID, IndX, IndY, GrowthDict, DispDict, GrowthList, MaintList, \
    DispList, ADList, EVList, TLList, RList, RVals, RX, RY, RID, RIDs, numDead = bide.reproduce(speciation, SpeciesIDs, Qs, IndIDs, IndID, \
    IndX, IndY,  width, height, GrowthDict, DispDict, RD, MaintDict, MainFactorDict, \
    RPFDict, EnvD, envgrads, GrowthList, MaintList, DispList, ADList, EVList, TLList, RList, RVals, RX, RY, RID, RIDs, \
    TrophicComplexityLevel, SpatialComplexityLevel, ResourceComplexityLevel, \
    BiologicalComplexityLevel, numDead)

    # maintenance
    SpeciesIDs, IndX, IndY, IndIDs, Qs, GrowthList, MaintList, DispList, ADList,\
    EVList, TLList, RList, RVals, RX, RY, RIDs, RID, numDead = bide.maintenance(SpeciesIDs, IndX, IndY, MaintDict, MainFactorDict, \
    RPFDict, EnvD, IndIDs, Qs, GrowthList, MaintList, DispList, ADList, EVList, TLList, RList, RVals, RX, RY, RIDs, RID,\
    TrophicComplexityLevel, SpatialComplexityLevel, ResourceComplexityLevel, \
    BiologicalComplexityLevel, numDead)

    # transition to or from dormancy
    SpeciesIDs, IndX, IndY, GrowthList, DispList, ADList, EVList, IndIDs, Qs, GrowthList, \
    MaintList, TLList, RList, RVals, RX, RY, RIDs, RID, numDead = bide.transition(SpeciesIDs, IndX, IndY, GrowthList, DispList, ADList, EVList,\
    IndIDs, Qs, MaintList, TLList, RList, RVals, RX, RY, RIDs, RID, MainFactorDict, RPFDict, \
    TrophicComplexityLevel, SpatialComplexityLevel, ResourceComplexityLevel, \
    BiologicalComplexityLevel, numDead)

    p2 = len(Qs)
    PRODI = p2 - p1

    ax = fig.add_subplot(111)
    plt.tick_params(axis='both', which='both', bottom='off', top='off', \
        left='off', right='off', labelbottom='off', labelleft='off')

    N = len(IndIDs)

    if N == 0 or N >= 1000:

        TrophicComplexityLevel = choice([1,2,3]) #
        SpatialComplexityLevel = choice([1,2,3]) # goes up to 3
        ResourceComplexityLevel = choice([1,2,3]) # goes up to 3 but needs enzyme breakdown
        BiologicalComplexityLevel = 2

        RDens, RDiv, RRich, Mu, Maint, ct, IndID, RID, N, T, R, PRODI, PRODQ = [0]*13
        ADList, ADs, AVG_DIST, SpecDisp, SpecMaint, SpecGrowth, GrowthList, MaintList, RVals, \
        DispList, EVList, TLList = [list([]) for _ in xrange(12)]

        SpeciesIDs, IndX, IndY, IndIDs, Qs, RX, RY, RIDs, RList, RVals, Gs, Ms, \
        Ds, Rs, PRODIs, Ns, RDENs, RDIVs, RRICHs, MUs, MAINTs, encList, Ragg, Iagg = [list([]) for _ in xrange(24)]

        p = 0
        BurnIn = 'not done'

    Ns.append(N)
    rr = len(RIDs)
    numD = ADList.count('d')

    pD = 0.0
    if N > 0:
        pD = round((numD/N*100), 2)

    Title = ['Active individuals (red) consume resources, grow, reproduce, go dormant (gray) and die in complex resource-limited environment.\n \
Resource complexity: ' + str(ResourceComplexityLevel) + ',  Trophic complexity: ' + str(TrophicComplexityLevel) + ',  \
Spatial complexity: ' + str(SpatialComplexityLevel) + '     N: '+str(N)+',  Resources: '+str(rr)+',  ct: '+str(ct)+ ', \
%Dormant: '+str(pD) + '\n# of inflowing resources: ' + str(r) + ', Aggregation: ' + str(round(std,2))]

    txt.set_text(' '.join(Title))
    ax.set_ylim(0, height)
    ax.set_xlim(0, width)

    if plot_system == 'yes':

        ##### PLOTTING THE SYSTEM ##############################################
        resource_scatImage.remove()
        Ind_scatImage.remove()
        colorlist = []

        ind_sizelist = []
        res_sizelist = []

        for val in RVals:
            res_sizelist.append(val*200)

        for i, val in enumerate(SpeciesIDs):
            if ADList[i] == 'a':
                colorlist.append('red')
            elif ADList[i] == 'd':
                colorlist.append('0.3')

            ind_sizelist.append(Qs[i] * 1000)


        resource_scatImage = ax.scatter(RX, RY, s = res_sizelist, c = 'w',
                    edgecolor = 'SpringGreen', lw = 2.0, alpha=0.7)

        Ind_scatImage = ax.scatter(IndX, IndY, s = ind_sizelist, c = colorlist,
                    edgecolor = '0.2', lw = 0.2, alpha=0.8)

    #if len(Ns) >= 100:
        #BurnIn = 'done'
    #    Ns = []

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

        PRODIs.append(PRODI)

        if len(RList) > 0:
            RDens = len(RList)/(height*width)

        RDENs.append(RDens)
        R = len(RX)
        Rs.append(R)
        encList.append(encounters)

        if N >= 1:

            if N >= 2:
                Imor = spatial.morisitas(IndX, IndY, width, height)
                Iagg.append(Imor)

            if R >= 1:
                q = min([20, R])

                #avg_dist1 = spatial.avg_dist(IndX, RX, IndY, RY, q)
                avg_dist2 = spatial.nearest_neighbor(IndX, RX, IndY, RY, q)
                AVG_DIST.append(avg_dist2)

                if R >= 2:
                    Rmor = spatial.morisitas(RX, RY, width, height)
                    Ragg.append(Rmor)

            spD = DispDict.values()
            spM = MaintDict.values()
            spG = GrowthDict.values()

            SpecDisp.append(mean(spD))
            SpecMaint.append(mean(spM))
            SpecGrowth.append(mean(spG))

            Gs.append(mean(GrowthList))
            Ms.append(mean(MaintList))
            Ds.append(mean(DispList))

            numD = ADList.count('d')
            ADs.append(numD/len(ADList))
            Deadlist.append(numDead)

        if len(Ns) >= 20:

            print sim, 'r:', r, 'R:', int(round(mean(Rs))), 'N:', int(round(mean(Ns))), \
            '%Dormant:', round(mean(ADs),3), 'Encounters:', round(mean(encList),2), 'Spatial:', SpatialComplexityLevel, \
            'Resource:', ResourceComplexityLevel, 'Trophic:', TrophicComplexityLevel, \
            'Agg(I):', round(mean(Iagg), 2), 'Agg(R):', round(mean(Ragg),2)

            OUT1 = open(GenPath + '/SimData.csv','a')

            outlist = [sim, mean(PRODIs), var(PRODIs), r, gmax, maintmax, dmax, alpha, seedCom, \
            width-0.2, height, mean(Ns), var(Ns), m, mean(RDENs), var(RDENs), mean(Rs), var(Rs), \
            speciation, mean(Gs), var(Gs), mean(Ms), var(Ms), mean(Ds), var(Ds), \
            mean(SpecGrowth), var(SpecGrowth), mean(SpecDisp), var(SpecDisp), \
            mean(SpecMaint),  var(SpecMaint), mean(AVG_DIST), var(AVG_DIST), \
            mean(ADs), var(ADs), TrophicComplexityLevel, ResourceComplexityLevel, \
            SpatialComplexityLevel, mean(encList), var(encList), std, mean(Iagg), \
            var(Iagg), mean(Ragg), var(Ragg), mean(Deadlist), var(Deadlist), ct]

            outlist = str(outlist).strip('[]')
            outlist = str(outlist).strip('')

            print>>OUT1, outlist
            OUT1.close()

            Rates = np.roll(Rates, -1, axis=0)
            u0 = Rates[0]

            if static == 'no':
                TrophicComplexityLevel = choice([1,2,3,4]) #
                SpatialComplexityLevel = choice([1,2,3]) # goes up to 3
                ResourceComplexityLevel = choice([1,2,3]) # goes up to 3 but needs enzyme breakdown
                BiologicalComplexityLevel = 2

            RDens, RDiv, RRich, Mu, Maint, ct, IndID, RID, N, T, R, PRODI, PRODQ, numD = [0]*14

            ADList, ADs, AVG_DIST, SpecDisp, SpecMaint, SpecGrowth, GrowthList, MaintList, RVals, \
            DispList, EVList, TLList, Deadlist = [list([]) for _ in xrange(13)]

            SpeciesIDs, IndX, IndY, IndIDs, Qs, RX, RY, RIDs, RList, RVals, Gs, Ms, \
            Ds, Rs, PRODIs, Ns, RDENs, RDIVs, RRICHs, MUs, MAINTs, encList, Ragg, Iagg = [list([]) for _ in xrange(24)]

            p = 0
            BurnIn = 'not done'

            if u0 == max(Rates):
                sim += 1

                width, height, alpha, speciation, seedCom, m, r, \
                gmax, maintmax, dmax, envgrads, Rates, pmax, mmax, std, maxgen = rp.get_rand_params(fixed)

                GrowthDict, MaintDict, MainFactorDict, RPFDict, EnvD, RD, DispDict,\
                EnvD = {}, {}, {}, {}, {}, {}, {}, {}

                enzyme_field = [0]*(width*height)

            ####################### REPLACE ENVIRONMENT ########################
            ax = fig.add_subplot(111)



############################# Conditions #######################################

fixed = True

TrophicComplexityLevel = choice([1,2,3,4]) #
SpatialComplexityLevel = choice([1,2,3]) # goes up to 3
ResourceComplexityLevel = choice([1,2,3]) # goes up to 3 but needs enzyme breakdown
BiologicalComplexityLevel = 2

################ Randomly chosen variables ##################################

width, height, alpha, speciation, seedCom, m, r, gmax, maintmax, dmax,\
envgrads, Rates, pmax, mmax, std, maxgen = rp.get_rand_params(fixed)

u0 = Rates[0]

#######################  Lists & Dictionaries  #########################
RDens, RDiv, RRich, Mu, Maint, ct, IndID, RID, N, T, R, PRODI, PRODQ, numD = [0]*14

ADList, ADs, AVG_DIST, SpecDisp, SpecMaint, SpecGrowth, GrowthList, MaintList, RList, \
DispList, EVList, TLList, Deadlist = [list([]) for _ in xrange(13)]

SpeciesIDs, IndX, IndY, IndIDs, Qs, RX, RY, RIDs, RList, RVals, Gs, Ms, \
Ds, Rs, PRODIs, Ns, RDENs, RDIVs, RRICHs, MUs, MAINTs, encList, Ragg, Iagg = [list([]) for _ in xrange(24)]

GrowthDict, MaintDict, MainFactorDict, RPFDict, EnvD, RD, DispDict, EnvD = {}, {}, {}, {}, {}, {}, {}, {}
enzyme_field = [0]*(width*height)

num_sims, LowerLimit, sim, p = 100, 30, 1, 0.0
static = 'no'
BurnIn = 'not done'

############### INITIALIZE GRAPHICS ############################################
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111) # initiate first plot
Ind_scatImage = ax.scatter([0],[0], alpha=0)
resource_scatImage = ax.scatter([0],[0], alpha=0)

Title = ['','']
txt = fig.suptitle(' '.join(Title), fontsize = 12)

ani = animation.FuncAnimation(fig, nextFrame, frames=400, interval=100, blit=False) # 20000 frames is a long movie
plt.show()
#ani.save(mydir+'/GitHub/Micro-Encounter/results/movies/Micro-Encounter.avi', bitrate=2000)
