from __future__ import division
import matplotlib.animation as animation
import matplotlib.pyplot as plt

import statsmodels.tsa.stattools as sta
from math import isnan

import numpy as np
from numpy import mean
import sys
import os

mydir = os.path.expanduser("~/")
sys.path.append(mydir + "GitHub/Micro-Encounter/model/metrics")
import metrics
sys.path.append(mydir + "GitHub/Micro-Encounter/model/bide")
import bide
sys.path.append(mydir + "GitHub/Micro-Encounter/model/randparams")
import randparams as rp
sys.path.append(mydir + "GitHub/Micro-Encounter/model/spatial")
import spatial
sys.path.append(mydir + "GitHub/Micro-Encounter/model/EnzymeField")
import EnzymeField as field


GenPath = mydir + 'GitHub/Micro-Encounter/results/simulated_data/'
OUT1 = open(GenPath + 'SimData.csv','w')

print>>OUT1, 'RowID, motion, ind.production, res.inflow, N.types, P.types, C.types, max.res.val, \
              max.growth.rate, max.met.maint, max.active.dispersal, logseries.a, \
              starting.seed, width, height, viscosity, total.abundance, immigration.rate, \
              resource.concentration, shannons.resource.diversity, resource.richness, \
              resource.particles, speciation, avg.per.capita.growth, avg.per.capita.maint, \
              avg.per.capita.N.efficiency, avg.per.capita.P.efficiency, \
              avg.per.capita.C.efficiency, avg.per.capita.active.dispersal, \
              amplitude, flux, frequency, phase, disturbance, spec.growth, \
              spec.disp, spec.maint, avg.dist, dorm.freq'
OUT1.close()


def nextFrame(arg):

    """ Function called for each successive animation frame; arg is the frame number """

    global mmax, pmax, ADs, ADList, AVG_DIST, SpecDisp, SpecMaint, SpecGrowth
    global fixed, p, BurnIn, t, num_sims, width, height, Rates, GrowthDict, N_RD
    global P_RD, C_RD, DispDict, MaintDict, gmax, dmax, maintmax, IndIDs, Qs
    global IndID, IndX, IndY, Ind_scatImage, SpeciesIDs
    global RTypes, RX, RY, RID, RIDs, RVals, EnvD, resource_scatImage, ct1, Mu, Maint
    global motion, speciation, seedCom, m, r, nNi, nP, nC, rmax, sim
    global N, ct, RDens, RDiv, RRich, T, R, LowerLimit, prod_i, prod_q, alpha
    global Ts, Rs, PRODIs, Ns, RDENs, RDIVs, RRICHs, GrowthList, MaintList, N_RList, P_RList
    global MUs, MAINTs, PRODNs, PRODPs, PRODCs, Gs, Ms, NRs, PRs, CRs, Ds
    global C_RList, DispList, envgrads, MainFactorDict, RPFDict, enzyme_field, u0

    ct += 1
    plot_system = 'no'

    # Inflow of resources
    RTypes, RVals, RX, RY,  RIDs, RID = bide.ResIn(RTypes, RVals, RX, RY,  RID, RIDs, r, rmax, nNi, nP, nC, width, height, u0)

    # Immigration
    SpeciesIDs, IndX, IndY,  MaintDict, MainFactorDict, RPFDict, EnvD, GrowthDict,\
    DispDict, IndIDs, IndID, Qs, N_RD, P_RD, C_RD, GrowthList, MaintList, N_RList, P_RList,\
    C_RList, DispList, ADList = bide.immigration(mmax, pmax, dmax, gmax, maintmax, seedCom, m, \
    SpeciesIDs, IndX, IndY, width, height, MaintDict, MainFactorDict, RPFDict, EnvD, envgrads,\
    GrowthDict, DispDict, IndIDs, IndID, Qs, N_RD, P_RD, C_RD, nNi, nP, nC, u0,\
    alpha, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList, ADList, ct)

    # dispersal
    #Lists = [SpeciesIDs, IndIDs, IndID, Qs, DispDict, GrowthList, MaintList, \
    #N_RList, P_RList, C_RList, DispList, ADList]

    #SpeciesIDs, IndX, IndY, IndIDs, IndID, Qs, GrowthList, \
    #MaintList, N_RList, P_RList, C_RList, DispList, ADList = bide.movement('individual',\
    #motion, Lists, IndX, IndY, width, height, u0)


    # Chemotaxis
    #SpeciesIDs, Qs, IndIDs, ID, X, Y, GrowthDict, DispDict, GrowthList, \
    #MaintList, N_RList, P_RList, C_RList, DispList, ADList = bide.chemotaxis(reproduction, \
    #speciation, SpeciesIDs, Qs, IndIDs, IndID, IndX, IndY, width, height, GrowthDict, \
    #DispDict, N_RD, P_RD, C_RD, MaintDict, MainFactorDict, RPFDict, EnvD, envgrads, nNi, \
    #nP, nC, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList, ADList)


    # Forage
    #SpeciesIDs, Qs, IndIDs, ID, X, Y, GrowthDict, DispDict, GrowthList, MaintList,\
    #N_RList, P_RList, C_RList, DispList, ADList = bide.density_forage(RVals, RX, RY, reproduction,\
    #speciation, SpeciesIDs, Qs, IndIDs, IndID, IndX, IndY,  width, height, GrowthDict,\
    #DispDict, N_RD, P_RD, C_RD, MaintDict, MainFactorDict, RPFDict, EnvD, envgrads, nNi,\
    #nP, nC, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList, ADList)

    PRODI = 0
    p1, TNQ1, TPQ1, TCQ1 = metrics.getprod(Qs)

    # Modify enzyme field
    enzyme_field = field.EnzymeField(enzyme_field, IndX, IndY, ADList, Qs, width)

    # Consume
    RTypes, RVals, RIDs, RID, RX, RY, SpeciesIDs, Qs, IndIDs, IndID, IndX, IndY, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList, ADList = bide.consume(enzyme_field, RTypes, RVals, RIDs, RID, RX, RY, SpeciesIDs, Qs, IndIDs, IndID, IndX, IndY,  width, height, GrowthDict, N_RD, P_RD, C_RD, DispDict, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList, ADList)

    # Reproduction
    SpeciesIDs, Qs, IndIDs, ID, X, Y, GrowthDict, DispDict, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList, ADList = bide.reproduce(speciation, SpeciesIDs, Qs, IndIDs, IndID, IndX, IndY,  width, height, GrowthDict, DispDict, N_RD, P_RD, C_RD, MaintDict, MainFactorDict, RPFDict, EnvD, envgrads, nNi, nP, nC, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList, ADList)

    # maintenance
    SpeciesIDs, X, Y, IndIDs, Qs, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList, ADList = bide.maintenance(SpeciesIDs, IndX, IndY, MaintDict, MainFactorDict, RPFDict, EnvD, IndIDs, Qs, GrowthList, MaintList, N_RList, P_RList, C_RList, DispList, ADList)

    # transition to or from dormancy
    Sp_IDs, IDs, Qs, GrowthList, MaintList, ADList = bide.transition(SpeciesIDs, IndIDs, Qs, GrowthList, MaintList, MainFactorDict, RPFDict,  ADList)

    p2, TNQ2, TPQ2, TCQ2 = metrics.getprod(Qs)
    PRODI = p2 - p1

    ax = fig.add_subplot(111)
    plt.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labelleft='off')

    N = len(IndIDs)
    rr = len(RIDs)
    numD = ADList.count('d')

    if N != len(ADList):
        print N, len(SpeciesIDs), len(ADList)
        print "N != len(ADList)"
        sys.exit()

    Title = ['Individuals consume resources, grow, reproduce, and die in complex resource-limited environment.\n \
    N: '+str(N)+', Resources: '+str(rr)+', ct: '+str(ct)+', %dormant: '+str(round((numD/N)*100, 2))]

    txt.set_text(' '.join(Title))
    ax.set_ylim(0, height)
    ax.set_xlim(0, width)

    if plot_system == 'yes':
        ##### PLOTTING THE SYSTEM ##############################################
        resource_scatImage.remove()
        Ind_scatImage.remove()
        colorlist = []
        sizelist = []
        for i, val in enumerate(SpeciesIDs):
            if ADList[i] == 'a':
                colorlist.append('red')
            elif ADList[i] == 'd':
                colorlist.append('0.3')

            sizelist.append(min(Qs[i]) * 1000)

        resource_scatImage = ax.scatter(RX, RY, s = RVals*100, c = 'w', edgecolor = 'SpringGreen', lw = 0.6, alpha=0.3)

        Ind_scatImage = ax.scatter(IndX, IndY, s = sizelist, c = colorlist, edgecolor = '0.2', lw = 0.2, alpha=0.9)

    Ns.append(N)

    if N == 0 and BurnIn == 'not done':
        Ns = [Ns[-1]] # only keep the most recent N value
        BurnIn = 'done'

    if len(Ns) >= 100:

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

        # Examining the resource RAD
        if len(RTypes) > 0:
            RRAD, Rlist = bide.GetRAD(RTypes)
            RDens = len(RTypes)/(height*width)
            RDiv = float(metrics.Shannons_H(RRAD))
            RRich = len(Rlist)

        RDENs.append(RDens)
        RDIVs.append(RDiv)
        RRICHs.append(RRich)

        R = len(RX)
        Rs.append(R)

        if N >= 1:

            if R >= 1:
                q = min([10, R])
                avg_dist = spatial.avg_dist(IndX, RX, IndY, RY, q)
                avg_dist = spatial.nearest_neighbor(IndX, RX, IndY, RY, q)

                AVG_DIST.append(avg_dist)

            spD = DispDict.values()
            spM = MaintDict.values()
            #spMF = MainFactorDict.values()
            #spRPF = RPFDict.values()
            spG = GrowthDict.values()

            SpecDisp.append(mean(spD))
            SpecMaint.append(mean(spM))
            SpecGrowth.append(mean(spG))

            Gs.append(mean(GrowthList))
            Ms.append(mean(MaintList))
            Ds.append(mean(DispList))

            numD = ADList.count('d')
            ADs.append(numD/len(ADList))


        if len(Ns) > 100:
            print sim, ' N:', int(round(mean(Ns))), 'dormant:', round(mean(ADs),3), 'distance:', round(mean(AVG_DIST))
            OUT1 = open(GenPath + '/SimData.csv','a')

            outlist = [sim, motion, mean(PRODIs), r, nNi, nP, nC, rmax, gmax, maintmax, dmax,\
            alpha, seedCom, u0, width-0.2, height, N, m, mean(RDENs), mean(RDIVs), mean(RRICHs),\
            T, R, speciation, mean(Gs), mean(Ms), mean(NRs), mean(PRs), mean(CRs), mean(Ds),\
            mean(SpecGrowth), mean(SpecDisp), mean(SpecMaint), mean(AVG_DIST), mean(ADs)]

            outlist = str(outlist).strip('[]')

            print>>OUT1, outlist
            OUT1.close()

            ct1 += 1
            ct = 0

            Rates = np.roll(Rates, -1, axis=0)
            u0 = Rates[0]

            RDens, RDiv, RRich, Mu, Maint, ct, IndID, RID, N, ct1, T, R, PRODI, PRODQ = [0]*14

            ADList, ADs, AVG_DIST, SpecDisp, SpecMaint, SpecGrowth, GrowthList, MaintList, N_RList, \
            P_RList, C_RList, DispList = [list([]) for _ in xrange(12)]

            SpeciesIDs, IndX, IndY, IndIDs, Qs, RX, RY, RIDs, RTypes, RVals, Gs, Ms, NRs, PRs, CRs,\
            Ds, Rs, PRODIs, Ns, RDENs, RDIVs, RRICHs, MUs, MAINTs = [list([]) for _ in xrange(24)]

            p = 0
            BurnIn = 'not done'

            if u0 == max(Rates):
                sim += 1

                width, height, alpha, motion, speciation, seedCom, m, r, nNi, nP, nC, rmax, gmax, maintmax, dmax, envgrads, Rates, pmax, mmax = rp.get_rand_params(fixed)
                GrowthDict, MaintDict, MainFactorDict, RPFDict, EnvD, N_RD, P_RD, C_RD, DispDict, EnvD = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}
                enzyme_field = [0]*(width*height)

            ####################### REPLACE ENVIRONMENT ########################
            ax = fig.add_subplot(111)


################ Randomly chosen variables ##################################
fixed = True
width, height, alpha, motion, speciation, seedCom, m, r, nNi, nP, nC, rmax, gmax, maintmax, dmax, envgrads, Rates, pmax, mmax = rp.get_rand_params(fixed)
u0 = Rates[0]

#######################  Lists & Dictionaries  #########################
RDens, RDiv, RRich, Mu, Maint, ct, IndID, RID, N, ct1, T, R, PRODI, PRODQ = [0]*14

ADList, ADs, AVG_DIST, SpecDisp, SpecMaint, SpecGrowth, GrowthList, MaintList, N_RList, \
P_RList, C_RList, DispList = [list([]) for _ in xrange(12)]

SpeciesIDs, IndX, IndY, IndIDs, Qs, RX, RY, RIDs, RTypes, RVals, Gs, Ms, NRs, PRs, CRs,\
Ds, Rs, PRODIs, Ns, RDENs, RDIVs, RRICHs, MUs, MAINTs = [list([]) for _ in xrange(24)]

GrowthDict, MaintDict, MainFactorDict, RPFDict, EnvD, N_RD, P_RD, C_RD, DispDict, EnvD = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}
enzyme_field = [0]*(width*height)

num_sims, LowerLimit, sim, p = 10000, 30, 1, 0.0
BurnIn = 'not done'

############### INITIALIZE GRAPHICS ############################################
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111) # initiate first plot
Ind_scatImage = ax.scatter([0],[0], alpha=0)
resource_scatImage = ax.scatter([0],[0], alpha=0)

Title = ['','']
txt = fig.suptitle(' '.join(Title), fontsize = 12)

ani = animation.FuncAnimation(fig, nextFrame, frames=110, interval=40, blit=False) # 20000 frames is a long movie
plt.show()
#ani.save(mydir+'/GitHub/Micro-Encounter/results/movies/examples/2015_10_05_1751_hydrobide.avi', bitrate=5000)
