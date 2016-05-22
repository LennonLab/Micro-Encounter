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
#dat = dat.convert_objects(convert_numeric=True).dropna()


#-------------------------DATA FILTERS------------------------------------------

#dat = dat[dat['SpatialComplexityLevel'] == 3]
#dat = dat[dat['ResourceComplexityLevel'] == 1]
dat = dat[dat['TrophicComplexityLevel'] != 3]
#dat = dat[dat['IncomingResAgg'] < 0.2]
#dat = dat[dat['MaxMetMaint'] < 0.005]
#dat = dat[dat['ResInflow'] < 18]

#-------------------------------------------------------------------------------

ComplexityLevels = ['res', 'troph', 'spatial']

for level in ComplexityLevels:

    #### plot figure ###############################################################
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

    ylab = 'Total abundance, '+r'$N$'
    xlab = 'Resource supply'

    plt.scatter(dat1['ResInflow'], dat1['MeanTotalAbundance'], color = 'm', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k', label=label1)
    plt.scatter(dat2['ResInflow'], dat2['MeanTotalAbundance'], color = 'steelblue', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k', label=label2)
    plt.scatter(dat3['ResInflow'], dat3['MeanTotalAbundance'], color = 'goldenrod', alpha = 0.8 , s = 10, linewidths = 0.0, edgecolor = 'k', label=label3)

    plt.ylabel(ylab, fontsize=fs+5)
    plt.xlabel(xlab, fontsize=fs+5)
    plt.ylim(1, 2000)
    plt.xlim(0.5, 100)
    plt.xscale('log')
    plt.yscale('log')
    plt.tick_params(axis='both', which='major', labelsize=fs)
    plt.legend(bbox_to_anchor=(-0.04, 1.05, 2.48, .2), loc=10, ncol=3, mode="expand",prop={'size':fs})


    #### PLOT 2 ################################
    fig.add_subplot(2, 2, 2)

    ylab = '% change in '+r'$N$'
    xlab = 'Resource supply'

    plt.scatter(dat1['ResInflow'], dat1['MeanIndProduction']/dat1['MeanTotalAbundance'], color = 'm', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
    plt.scatter(dat2['ResInflow'], dat2['MeanIndProduction']/dat2['MeanTotalAbundance'], color = 'steelblue', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
    plt.scatter(dat3['ResInflow'], dat3['MeanIndProduction']/dat3['MeanTotalAbundance'], color = 'goldenrod', alpha = 0.8 , s = 10, linewidths = 0.0, edgecolor = 'k')

    plt.ylabel(ylab, fontsize=fs+5)
    plt.xlabel(xlab, fontsize=fs+5)
    plt.xscale('log')
    plt.ylim(-0.0, 0.1)
    plt.xlim(0.5, 150)
    plt.tick_params(axis='both', which='major', labelsize=fs)

    #### PLOT 3 #################################################################
    fig.add_subplot(2, 2, 3)

    ylab = 'Variance in '+r'$N$'
    xlab = 'Resource supply'

    plt.scatter(dat1['ResInflow'], dat1['VarTotalAbundance'], color = 'm', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
    plt.scatter(dat2['ResInflow'], dat2['VarTotalAbundance'], color = 'steelblue', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
    plt.scatter(dat3['ResInflow'], dat3['VarTotalAbundance'], color = 'goldenrod', alpha = 0.8 , s = 10, linewidths = 0.0, edgecolor = 'k')

    plt.ylabel(ylab, fontsize=fs+5)
    plt.xlabel(xlab, fontsize=fs+5)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(0.1, 20000)
    plt.xlim(0.5, 150)
    plt.tick_params(axis='both', which='major', labelsize=fs)


    #### PLOT 4 #################################################################
    fig.add_subplot(2, 2, 4)

    ylab = 'Individual maintenance'
    xlab = 'Resource supply'

    plt.scatter(dat1['ResInflow'], dat1['MeanPerCapitaMaint'], color = 'm', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
    plt.scatter(dat2['ResInflow'], dat2['MeanPerCapitaMaint'], color = 'steelblue', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
    plt.scatter(dat3['ResInflow'], dat3['MeanPerCapitaMaint'], color = 'goldenrod', alpha = 0.8 , s = 10, linewidths = 0.0, edgecolor = 'k')

    plt.ylabel(ylab, fontsize=fs+5)
    plt.xlabel(xlab, fontsize=fs+5)
    plt.xscale('log')
    plt.ylim(0.0001, 0.0013)
    plt.xlim(0.5, 150)
    plt.tick_params(axis='both', which='major', labelsize=fs)

    #### Final Format and Save #####################################################
    plt.subplots_adjust(wspace=0.4, hspace=0.4)

    if level == 'spatial':
        plt.savefig(mydir + '/results/figures/ResInflow-SpatialComplexity-NoTrophicLevel3.png', dpi=600, bbox_inches = "tight")
    elif level == 'res':
        plt.savefig(mydir + '/results/figures/ResInflow-ResourceComplexity-NoTrophicLevel3.png', dpi=600, bbox_inches = "tight")
    elif level == 'troph':
        plt.savefig(mydir + '/results/figures/ResInflow-TrophicComplexity-NoTrophicLevel3.png', dpi=600, bbox_inches = "tight")

    #plt.show()
