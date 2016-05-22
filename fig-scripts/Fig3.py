from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys

#import statsmodels.stats.api as sms
#import statsmodels.api as sm
import statsmodels.formula.api as smf
#from statsmodels.sandbox.regression.predstd import wls_prediction_std
from statsmodels.stats.outliers_influence import summary_table

complexity = 'res'
complexity = 'troph'
complexity = 'spatial'

mydir = os.path.expanduser('~/GitHub/Micro-Encounter')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")
dat = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')


#### plot figure ###############################################################
fs = 8 # fontsize
fig = plt.figure()


if complexity == 'res':

    dat1 = dat[dat['ResourceComplexityLevel'] == 1]
    dat2 = dat[dat['ResourceComplexityLevel'] == 2]
    dat3 = dat[dat['ResourceComplexityLevel'] == 3]

    label1 = 'No diversity, No complexity'
    label2 = 'No complexity'
    label3 = 'Diversity + Complexity'

if complexity == 'troph':

    dat1 = dat[dat['TrophicComplexityLevel'] == 1]
    dat2 = dat[dat['TrophicComplexityLevel'] == 2]

    label1 = 'No trophic complexity'
    label2 = 'Trophic complexity'

if complexity == 'spatial':

    dat1 = dat[dat['SpatialComplexityLevel'] == 1]
    dat2 = dat[dat['SpatialComplexityLevel'] == 2]
    dat3 = dat[dat['SpatialComplexityLevel'] == 3]

    label1 = 'White noise'
    label2 = 'Aggregated w/ Random walks'
    label3 = 'Aggregated w/ chemotaxis'


#### PLOT 1 #################################################################
fig.add_subplot(2, 2, 1)

ylab = 'Percent Dormant'
xlab = 'Active dispersal'

plt.scatter(dat1['MeanPerCapitaActiveDispersal'], dat1['MeanDormFreq'], color = 'm', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k', label=label1)
plt.scatter(dat2['MeanPerCapitaActiveDispersal'], dat2['MeanDormFreq'], color = 'steelblue', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k', label=label2)
if complexity != 'troph':
    plt.scatter(dat3['MeanPerCapitaActiveDispersal'], dat3['MeanDormFreq'], color = 'goldenrod', alpha = 0.8 , s = 10, linewidths = 0.0, edgecolor = 'k', label=label3)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
#plt.yscale('log')
plt.xlim(0.0, 0.001)
plt.tick_params(axis='both', which='major', labelsize=fs)
plt.legend(bbox_to_anchor=(-0.04, 1.05, 2.48, .2), loc=10, ncol=3, mode="expand",prop={'size':fs})


#### PLOT 2 ################################
fig.add_subplot(2, 2, 2)

ylab = 'Mean Encounters'
xlab = 'Active dispersal'

plt.scatter(dat1['MeanPerCapitaActiveDispersal'], dat1['MeanEncounter'], color = 'm', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
plt.scatter(dat2['MeanPerCapitaActiveDispersal'], dat2['MeanEncounter'], color = 'steelblue', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
if complexity != 'troph':
    plt.scatter(dat3['MeanPerCapitaActiveDispersal'], dat3['MeanEncounter'], color = 'goldenrod', alpha = 0.8 , s = 10, linewidths = 0.0, edgecolor = 'k')

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
#plt.xscale('log')
plt.yscale('log')
plt.ylim(0.008, 1000)
plt.xlim(0.0, 0.001)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### PLOT 3 #################################################################
fig.add_subplot(2, 2, 3)

ylab = 'Mean Distance'
xlab = 'Active dispersal'

plt.scatter(dat1['MeanPerCapitaActiveDispersal'], dat1['MeanAvgDist'], color = 'm', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
plt.scatter(dat2['MeanPerCapitaActiveDispersal'], dat2['MeanAvgDist'], color = 'steelblue', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
if complexity != 'troph':
    plt.scatter(dat3['MeanPerCapitaActiveDispersal'], dat3['MeanAvgDist'], color = 'goldenrod', alpha = 0.8 , s = 10, linewidths = 0.0, edgecolor = 'k')

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
#plt.yscale('log')
#plt.xlim(1, 300)
plt.xlim(0.0, 0.001)
plt.tick_params(axis='both', which='major', labelsize=fs)


#### PLOT 4 #################################################################
fig.add_subplot(2, 2, 4)

ylab = 'Total abundance'
xlab = 'Active dispersal'

plt.scatter(dat1['MeanPerCapitaActiveDispersal'], dat1['MeanTotalAbundance'], color = 'm', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
plt.scatter(dat2['MeanPerCapitaActiveDispersal'], dat2['MeanTotalAbundance'], color = 'steelblue', alpha = 0.7 , s = 10, linewidths = 0.0, edgecolor = 'k')
if complexity != 'troph':
    plt.scatter(dat3['MeanPerCapitaActiveDispersal'], dat3['MeanTotalAbundance'], color = 'goldenrod', alpha = 0.8 , s = 10, linewidths = 0.0, edgecolor = 'k')

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.yscale('log')
plt.xlim(0.0, 0.001)
#plt.ylim(1, 2000)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)

if complexity == 'spatial':
    plt.savefig(mydir + '/results/figures/Fig3-SpatialComplexity.png', dpi=600, bbox_inches = "tight")
elif complexity == 'res':
    plt.savefig(mydir + '/results/figures/Fig3-ResourceComplexity.png', dpi=600, bbox_inches = "tight")
elif complexity == 'troph':
    plt.savefig(mydir + '/results/figures/Fig3-TrophicComplexity.png', dpi=600, bbox_inches = "tight")

#plt.show()
