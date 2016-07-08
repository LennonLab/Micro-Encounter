from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys

mydir = os.path.expanduser('~/GitHub/Micro-Encounter')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")
dat = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')
dat = dat.convert_objects(convert_numeric=True).dropna()

#print dat.shape
#sys.exit()

#-------------------------DATA FILTERS------------------------------------------

dat = dat[dat['ResourceComplexityLevel'] == 1]
dat = dat[dat['TrophicComplexityLevel'] == 1]

#-------------------------------------------------------------------------------

ComplexityLevels = ['spatial']

for level in ComplexityLevels:


    #### plot figure ###############################################################
    fs = 8 # fontsize
    fig = plt.figure()

    if level == 'spatial':

        dat1 = dat[dat['SpatialComplexityLevel'] == 1]
        dat2 = dat[dat['SpatialComplexityLevel'] == 2]
        dat3 = dat[dat['SpatialComplexityLevel'] == 3]

        label1 = 'White noise'
        label2 = 'Aggregated w/ Random walks'
        label3 = 'Aggregated w/ chemotaxis'


    #### PLOT 1 #################################################################
    fig.add_subplot(2, 2, 1)

    xlab = 'Average encounters'
    ylab = '% Dormancy'

    plt.scatter(dat1['MeanEncounter'], dat1['MeanDormFreq']*100, color = 'm', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k', label=label1)
    plt.scatter(dat2['MeanEncounter'], dat2['MeanDormFreq']*100, color = 'steelblue', alpha = 0.8 , s = 10, linewidths = 0.0, edgecolor = 'k', label=label2)
    plt.scatter(dat3['MeanEncounter'], dat3['MeanDormFreq']*100, color = 'goldenrod', alpha = 0.8 , s = 10, linewidths = 0.0, edgecolor = 'k', label=label3)

    plt.ylabel(ylab, fontsize=fs+5)
    plt.xlabel(xlab, fontsize=fs+5)
    plt.xscale('log')
    plt.xlim(0.1, 300)
    plt.tick_params(axis='both', which='major', labelsize=fs)
    plt.legend(bbox_to_anchor=(-0.04, 1.05, 2.48, .2), loc=10, ncol=3, mode="expand",prop={'size':fs})


    #### PLOT 2 ################################
    fig.add_subplot(2, 2, 2)

    xlab = 'Average encounters'
    ylab = 'Individual production'

    plt.scatter(dat1['MeanEncounter'], dat1['MeanIndProduction'], color = 'm', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
    plt.scatter(dat2['MeanEncounter'], dat2['MeanIndProduction'], color = 'steelblue', alpha = 0.8 , s = 10, linewidths = 0.0, edgecolor = 'k')
    plt.scatter(dat3['MeanEncounter'], dat3['MeanIndProduction'], color = 'goldenrod', alpha = 0.8 , s = 10, linewidths = 0.0, edgecolor = 'k')

    plt.ylabel(ylab, fontsize=fs+5)
    plt.xlabel(xlab, fontsize=fs+5)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(0.008, 90)
    plt.xlim(0.1, 300)
    plt.tick_params(axis='both', which='major', labelsize=fs)

    #### PLOT 3 #################################################################
    fig.add_subplot(2, 2, 3)

    xlab = 'Average encounters'
    ylab = 'Resource aggregation'

    plt.scatter(dat1['MeanEncounter'], dat1['MeanTotalAbundance'], color = 'm', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
    plt.scatter(dat2['MeanEncounter'], dat2['MeanTotalAbundance'], color = 'steelblue', alpha = 0.8 , s = 10, linewidths = 0.0, edgecolor = 'k')
    plt.scatter(dat3['MeanEncounter'], dat3['MeanTotalAbundance'], color = 'goldenrod', alpha = 0.8 , s = 10, linewidths = 0.0, edgecolor = 'k')

    plt.ylabel(ylab, fontsize=fs+5)
    plt.xlabel(xlab, fontsize=fs+5)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(0.15, 300)
    plt.ylim(1, 3000)
    plt.tick_params(axis='both', which='major', labelsize=fs)


    #### PLOT 4 #################################################################
    fig.add_subplot(2, 2, 4)

    xlab = 'Average encounters'
    ylab = 'Active abundance'

    plt.scatter(dat1['MeanEncounter'], dat1['MeanTotalAbundance'] * (1 - dat1['MeanDormFreq']), color = 'm', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
    plt.scatter(dat2['MeanEncounter'], dat2['MeanTotalAbundance'] * (1 - dat2['MeanDormFreq']), color = 'steelblue', alpha = 0.8 , s = 10, linewidths = 0.0, edgecolor = 'k')
    plt.scatter(dat3['MeanEncounter'], dat3['MeanTotalAbundance'] * (1 - dat3['MeanDormFreq']), color = 'goldenrod', alpha = 0.8 , s = 10, linewidths = 0.0, edgecolor = 'k')

    plt.ylabel(ylab, fontsize=fs+5)
    plt.xlabel(xlab, fontsize=fs+5)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(0.15, 1000)
    plt.ylim(0.9, 1000)
    plt.tick_params(axis='both', which='major', labelsize=fs)

    #### Final Format and Save #####################################################
    plt.subplots_adjust(wspace=0.4, hspace=0.4)

    if level == 'spatial':
        plt.savefig(mydir + '/results/figures/Encounters-SpatialComplexity.png', dpi=600, bbox_inches = "tight")
    elif level == 'res':
        plt.savefig(mydir + '/results/figures/Encounters-ResourceComplexity.png', dpi=600, bbox_inches = "tight")
    elif level == 'troph':
        plt.savefig(mydir + '/results/figures/Encounters-TrophicComplexity.png', dpi=600, bbox_inches = "tight")

    #plt.show()
    plt.close()
