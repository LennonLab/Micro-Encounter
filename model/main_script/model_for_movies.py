from __future__ import division
import matplotlib.animation as animation
import matplotlib.pyplot as plt

import statsmodels.tsa.stattools as sta
from random import choice
from math import isnan

import numpy as np
from numpy import mean
import sys
import os

mydir = os.path.expanduser("~/")
sys.path.append(mydir + "GitHub/Micro-Encounter/model/bide")
import bide
sys.path.append(mydir + "GitHub/Micro-Encounter/model/randparams")
import randparams as rp
sys.path.append(mydir + "GitHub/Micro-Encounter/model/spatial")
import spatial


def nextFrame(arg):

    """ Function called for each successive animation frame; arg is the frame number """

    global mmax, pmax, ADs, ADList, AVG_DIST, SpecDisp, SpecMaint, SpecGrowth
    global fixed, p, BurnIn, t, num_sims, width, height, Rates, GrowthDict, RD
    global DispDict, MaintDict, gmax, dmax, maintmax, IndIDs, Qs, EVList
    global IndID, IndX, IndY, Ind_scatImage, SpeciesIDs, TLList
    global RX, RY, RID, RIDs, RVals, EnvD, resource_scatImage, Mu, Maint
    global seedCom, m, r, sim
    global N, ct, RDens, RDiv, RRich, T, R, LowerLimit, prod_i, prod_q
    global Ts, Rs, PRODIs, Ns, RDENs, RDIVs, RRICHs, GrowthList, MaintList, RList
    global MUs, MAINTs, PRODNs, PRODPs, PRODCs, Gs, Ms, Ds
    global DispList, envgrads, MainFactorDict, RPFDict, enzyme_field, u0

    global TrophicComplexityLevel, SpatialComplexityLevel, encList, std
    global ResourceComplexityLevel, BiologicalComplexityLevel
    global Ragg, Iagg, static, Deadlist

    numDead = 0
    ct += 1
    #print ct

    encounters = 0
    # Inflow of resources
    RList, RVals, RX, RY,  RIDs, RID = bide.ResIn(std, ct, RList, RVals, RX, RY,  RID,\
    RIDs, r, width, height, u0, TrophicComplexityLevel, SpatialComplexityLevel, \
    ResourceComplexityLevel, BiologicalComplexityLevel)

    # Immigration
    SpeciesIDs, IndX, IndY,  MaintDict, MainFactorDict, RPFDict, EnvD, GrowthDict,\
    DispDict, IndIDs, IndID, Qs, RD, GrowthList, MaintList, DispList, ADList, EVList,\
    TLList = bide.immigration(std, mmax, pmax, dmax, gmax, maintmax, seedCom, m, \
    SpeciesIDs, IndX, IndY, width, height, MaintDict, MainFactorDict, RPFDict, EnvD, envgrads,\
    GrowthDict, DispDict, IndIDs, IndID, Qs, RD, u0, GrowthList, MaintList, \
    DispList, ADList, EVList, TLList, ct, TrophicComplexityLevel, \
    SpatialComplexityLevel, ResourceComplexityLevel, BiologicalComplexityLevel)

    # Dispersal
    if SpatialComplexityLevel < 3:
        SpeciesIDs, Qs, IndIDs, ID, IndX, IndY, GrowthDict, DispDict, GrowthList, \
        MaintList, DispList, ADList, EVList, TLList, RList, RVals, RX, RY, RIDs, RID, numDead = bide.dispersal(SpeciesIDs,\
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

    # Consume
    RList, RVals, RIDs, RID, RX, RY, SpeciesIDs, Qs, encounters = bide.consume(enzyme_field, \
    RList, RVals, RIDs, RID, RX, RY, SpeciesIDs, Qs, IndIDs, IndID, IndX, IndY, \
    width, height, GrowthDict, MaintDict, RD, DispDict, GrowthList, MaintList, DispList, ADList, \
    TLList, EVList, TrophicComplexityLevel, SpatialComplexityLevel, \
    ResourceComplexityLevel, BiologicalComplexityLevel)

    # Reproduction
    SpeciesIDs, Qs, IndIDs, ID, IndX, IndY, GrowthDict, DispDict, GrowthList, MaintList, \
    DispList, ADList, EVList, TLList, RList, RVals, RX, RY, RID, RIDs, numDead = bide.reproduce(SpeciesIDs, Qs, IndIDs, IndID, \
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

    ax1 = fig.add_subplot(2,2,1) # initiate 1st plot
    ax2 = fig.add_subplot(2,2,2) # initiate 2nd plot
    ax3 = fig.add_subplot(2,2,3) # initiate 3rd plot
    ax4 = fig.add_subplot(2,2,4) # initiate 4th plot

    plt.tick_params(axis='both', which='both', bottom='off', top='off', \
        left='off', right='off', labelbottom='off', labelleft='off')

    N = len(IndIDs)

    if N == 0 or N > 20000:

        TrophicComplexityLevel = choice([1,2,3]) #
        SpatialComplexityLevel = choice([1,2,3]) # goes up to 3
        ResourceComplexityLevel = choice([1,2,3]) # goes up to 3 but needs enzyme breakdown
        BiologicalComplexityLevel = 2

        RDens, RDiv, RRich, Mu, Maint, ct, IndID, RID, N, T, R, PRODI, PRODQ = [0]*13
        ADList, ADs, AVG_DIST, SpecDisp, SpecMaint, SpecGrowth, GrowthList, MaintList, RVals, \
        DispList, EVList, TLList = [list([]) for _ in xrange(12)]

        SpeciesIDs, IndX, IndY, IndIDs, Qs, RX, RY, RIDs, RList, RVals, Gs, Ms, \
        Ds, Rs, PRODIs, Ns, RDENs, RDIVs, RRICHs, MUs, MAINTs, encList, Ragg, Iagg = [list([]) for _ in xrange(24)]
        
        if u0 == max(Rates):
            #sim += 1

            width, height, seedCom, m, r, gmax, maintmax, dmax, envgrads, Rates, pmax, mmax, std = rp.get_rand_params(fixed)

            GrowthDict, MaintDict, MainFactorDict, RPFDict, EnvD, RD, DispDict,\
            EnvD = {}, {}, {}, {}, {}, {}, {}, {}

            enzyme_field = [0]*(width*height)

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
    #ax1.set_ylim(0, height)
    #ax1.set_xlim(0, width)

    plot_system = 'yes'
    if plot_system == 'yes':

        ##### PLOTTING THE SYSTEM ##############################################
        resource_scatImage.remove()
        Ind_scatImage.remove()
        plot2.remove()
        plot3.remove()
        plot4.remove()
        
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


        resource_scatImage = ax1.scatter(RX, RY, s = res_sizelist, c = 'w',
                    edgecolor = 'SpringGreen', lw = 2.0, alpha=0.7)

        Ind_scatImage = ax1.scatter(IndX, IndY, s = ind_sizelist, c = colorlist,
                    edgecolor = '0.2', lw = 0.2, alpha=0.8)
                    
        ax2 = fig.add_subplot(2,2,2) # initiate 2nd plot
        ax3 = fig.add_subplot(2,2,3) # initiate 3rd plot
        ax4 = fig.add_subplot(2,2,4) # initiate 4th plot

        plot2 = ax2.scatter([0],[0], alpha=0)
        plot3 = ax3.scatter([0],[0], alpha=0)
        plot4 = ax4.scatter([0],[0], alpha=0)

    if len(Ns) >= 50:

        if BurnIn == 'not done':
            AugmentedDickeyFuller = sta.adfuller(Ns)
            val, p = AugmentedDickeyFuller[0:2]

            if p >= 0.05:
                Ns.pop(0)

            elif p < 0.05 or isnan(p) == True:
                BurnIn = 'done'
                Ns = [Ns[-1]] # only keep the most recent N value
                
    if ct == 100:
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

            print '%4s' % sim, ' r:','%4s' %  r, ' R:','%4s' % int(round(mean(Rs))), ' N:','%5s' % int(round(mean(Ns))), \
            ' Dormant:', '%5s' % round(mean(ADs),3), ' Encounters:','%5s' % round(mean(encList),2), '   Spatial:', SpatialComplexityLevel, \
            'Resource:', ResourceComplexityLevel, 'Trophic:', TrophicComplexityLevel#, \
            #'Agg(I):', round(mean(Iagg), 2), 'Agg(R):', round(mean(Ragg),2)

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

                width, height, seedCom, m, r, gmax, maintmax, dmax, envgrads, Rates, pmax, mmax, std = rp.get_rand_params(fixed)

                GrowthDict, MaintDict, MainFactorDict, RPFDict, EnvD, RD, DispDict,\
                EnvD = {}, {}, {}, {}, {}, {}, {}, {}

                enzyme_field = [0]*(width*height)

            ####################### REPLACE ENVIRONMENT ########################
            ax1 = fig.add_subplot(2,2,1)
            ax2 = fig.add_subplot(2,2,2) # initiate first plot
            ax3 = fig.add_subplot(2,2,3) # initiate first plot
            ax4 = fig.add_subplot(2,2,4) # initiate first plot
            
            
############################# Conditions #######################################

fixed = True

TrophicComplexityLevel = choice([1,2,3,4]) #
SpatialComplexityLevel = choice([1,2,3]) # goes up to 3
ResourceComplexityLevel = choice([1,2,3]) # goes up to 3 but needs enzyme breakdown
BiologicalComplexityLevel = 2

################ Randomly chosen variables ##################################

width, height, seedCom, m, r, gmax, maintmax, dmax, envgrads, Rates, pmax, mmax, std = rp.get_rand_params(fixed)

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
fig = plt.figure()

ax1 = fig.add_subplot(2,2,1) # initiate 1st plot
ax2 = fig.add_subplot(2,2,2) # initiate 2nd plot
ax3 = fig.add_subplot(2,2,3) # initiate 3rd plot
ax4 = fig.add_subplot(2,2,4) # initiate 4th plot

Ind_scatImage = ax1.scatter([0],[0], alpha=0)
resource_scatImage = ax1.scatter([0],[0], alpha=0)
plot2 = ax2.scatter([0],[0], alpha=0)
plot3 = ax3.scatter([0],[0], alpha=0)
plot4 = ax4.scatter([0],[0], alpha=0)

Title = ['','']
txt = fig.suptitle(' '.join(Title), fontsize = 12)

ani = animation.FuncAnimation(fig, nextFrame, frames=400, interval=100, blit=False) # 20000 frames is a long movie
plt.show()
#ani.save(mydir+'/GitHub/Micro-Encounter/results/movies/Micro-Encounter.avi', bitrate=2000)
