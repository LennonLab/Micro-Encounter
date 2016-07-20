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


#-------------------------DATA FILTERS------------------------------------------

#dat = dat[dat['SpatialComplexityLevel'] == 3]
#dat = dat[dat['ResourceComplexityLevel'] == 1]
dat = dat[dat['TrophicComplexityLevel'] != 2]
#dat = dat[dat['IncomingResAgg'] < 0.2]
#dat = dat[dat['MaxMetMaint'] < 0.005]
#dat = dat[dat['ResInflow'] < 18]

#-------------------------------------------------------------------------------

ComplexityLevels = ['res', 'troph', 'spatial']

for level in ComplexityLevels:


    #### plot figure ###########################################################

    lims = 'y'
    fs = 8 # fontsize
    fig = plt.figure()

    if level == 'res':

        dat1 = dat[dat['ResourceComplexityLevel'] == 1]
        dat2 = dat[dat['ResourceComplexityLevel'] == 2]
        dat3 = dat[dat['ResourceComplexityLevel'] == 3]

        label1 = 'No diversity, No complexity'
        label2 = 'No complexity'
        label3 = 'Diversity + Complexity'

    if level == 'troph':

        dat1 = dat[dat['TrophicComplexityLevel'] == 1]
        dat2 = dat[dat['TrophicComplexityLevel'] == 2]
        dat3 = dat[dat['TrophicComplexityLevel'] == 3]

        label1 = 'No trophic complexity'
        label2 = 'Trophic complexity'
        label3 = 'Recycling'

    if level == 'spatial':

        dat1 = dat[dat['SpatialComplexityLevel'] == 1]
        dat2 = dat[dat['SpatialComplexityLevel'] == 2]
        dat3 = dat[dat['SpatialComplexityLevel'] == 3]

        label1 = 'White noise'
        label2 = 'Aggregated w/ Random walks'
        label3 = 'Aggregated w/ chemotaxis'


    #### PLOT 1 #################################################################
    fig.add_subplot(2, 2, 1)

    xlab = 'Resource concentration'
    ylab = 'Resource aggregation'

    plt.scatter(dat1['MeanResourceConcentration'], dat1['MeanResAgg'], color = 'm', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k', label=label1)
    plt.scatter(dat2['MeanResourceConcentration'], dat2['MeanResAgg'], color = 'steelblue', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k', label=label2)
    plt.scatter(dat3['MeanResourceConcentration'], dat3['MeanResAgg'], color = 'goldenrod', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k', label=label3)

    plt.ylabel(ylab, fontsize=fs+5)
    plt.xlabel(xlab, fontsize=fs+5)

    plt.yscale('log')
    plt.xscale('log')

    if lims == 'y':
        plt.ylim(0.1, 1000)
        plt.xlim(0.001, 4.0)
    plt.tick_params(axis='both', which='major', labelsize=fs)
    plt.legend(bbox_to_anchor=(-0.04, 1.05, 2.48, .2), loc=10, ncol=3, mode="expand",prop={'size':fs})

    '''
    a = plt.axes([0.17, 0.78, 0.1, 0.1], axisbg='w')
    plt.scatter(dat1['VarResourceConcentration'], dat1['VarIndAgg'], color = 'm', alpha = 0.7 )
    plt.scatter(dat2['VarResourceConcentration'], dat2['VarIndAgg'], color = 'steelblue', alpha = 0.8)
    plt.scatter(dat3['VarResourceConcentration'], dat3['VarIndAgg'], color = 'goldenrod', alpha = 0.8)
    #plt.set_xlim(0.5*min(dat1['VarResourceConcentration']), 2*max(dat1['VarResourceConcentration']))
    #plt.set_ylim(0.5*min(dat1['VarIndAgg']), 2*max(dat1['VarIndAgg']))
    plt.title('Probability')
    plt.xticks([])
    plt.yticks([])
    '''

    #### PLOT 2 ################################
    fig.add_subplot(2, 2, 2)

    ylab = 'Avg Encounters'
    xlab = 'Resource aggregation'

    plt.scatter(dat1['MeanResAgg'], dat1['MeanEncounter'], color = 'm', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
    plt.scatter(dat2['MeanResAgg'], dat2['MeanEncounter'], color = 'steelblue', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
    plt.scatter(dat3['MeanResAgg'], dat3['MeanEncounter'], color = 'goldenrod', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')

    plt.ylabel(ylab, fontsize=fs+5)
    plt.xlabel(xlab, fontsize=fs+5)

    plt.xscale('log')
    plt.yscale('log')

    if lims == 'y':
        plt.ylim(0.01, 60)
        plt.xlim(0.1, 1000)

    plt.tick_params(axis='both', which='major', labelsize=fs)


    #### PLOT 3 #################################################################
    fig.add_subplot(2, 2, 3)

    ylab = 'Percent Dormant'
    xlab = 'Resource aggregation'

    plt.scatter(dat1['MeanResAgg'], dat1['MeanDormFreq'], color = 'm', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
    plt.scatter(dat2['MeanResAgg'], dat2['MeanDormFreq'], color = 'steelblue', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
    plt.scatter(dat3['MeanResAgg'], dat3['MeanDormFreq'], color = 'goldenrod', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')

    plt.ylabel(ylab, fontsize=fs+5)
    plt.xlabel(xlab, fontsize=fs+5)

    plt.xscale('log')

    if lims == 'y':
        plt.ylim(0.0, 1.05)
        plt.xlim(0.1, 1000)

    plt.tick_params(axis='both', which='major', labelsize=fs)


    #### PLOT 4 #################################################################
    fig.add_subplot(2, 2, 4)

    ylab = 'Productivity'
    xlab = 'Resource aggregation'

    plt.scatter(dat1['MeanResAgg'], dat1['MeanIndProduction'], color = 'm', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
    plt.scatter(dat2['MeanResAgg'], dat2['MeanIndProduction'], color = 'steelblue', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
    plt.scatter(dat3['MeanResAgg'], dat3['MeanIndProduction'], color = 'goldenrod', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')

    plt.ylabel(ylab, fontsize=fs+5)
    plt.xlabel(xlab, fontsize=fs+5)

    plt.yscale('log')
    plt.xscale('log')

    if lims == 'y':
        plt.ylim(0.01, 40)
        plt.xlim(0.1, 1000)
    plt.tick_params(axis='both', which='major', labelsize=fs)

    #### Final Format and Save #####################################################
    plt.subplots_adjust(wspace=0.4, hspace=0.4)

    if level == 'spatial':
        plt.savefig(mydir + '/results/figures/Aggregation-SpatialComplexity-NoTrophicLevel2.png', dpi=600, bbox_inches = "tight")
    elif level == 'res':
        plt.savefig(mydir + '/results/figures/Aggregation-ResourceComplexity-NoTrophicLevel2.png', dpi=600, bbox_inches = "tight")
    elif level == 'troph':
        plt.savefig(mydir + '/results/figures/Aggregation-TrophicComplexity-NoTrophicLevel2.png', dpi=600, bbox_inches = "tight")

    #plt.show()
    plt.close()
