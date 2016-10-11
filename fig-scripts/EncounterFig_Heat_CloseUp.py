from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys

mydir = os.path.expanduser('~/GitHub/Micro-Encounter')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")
df = pd.read_csv(mydir + '/results/simulated_data/2016_09_17_SimData.csv')

#-------------------------DATA TRANSFORMATIONS -----------------------
df2 = pd.DataFrame({'Encounters' : np.log10(df['Encounters'].groupby(df['sim']).mean())})

df2['SpatialComplexity'] = df['SpatialComplexity'].groupby(df['sim']).unique()
df2['SpatialComplexity'] = [df2['SpatialComplexity'][i][0] for i in df2['SpatialComplexity'].keys()]

df2['TrophicComplexity'] = df['TrophicComplexity'].groupby(df['sim']).unique()
df2['TrophicComplexity'] = [df2['TrophicComplexity'][i][0] for i in df2['TrophicComplexity'].keys()]

df2['ResourceComplexity'] = df['ResourceComplexity'].groupby(df['sim']).unique()
df2['ResourceComplexity'] = [df2['ResourceComplexity'][i][0] for i in df2['ResourceComplexity'].keys()]

df2['D'] = df['numDead'].groupby(df['sim']).mean()
df2['R'] = df['R'].groupby(df['sim']).mean()
df2['Q'] = df['MeanCellQuota'].groupby(df['sim']).mean()
df2['N'] = df['N'].groupby(df['sim']).mean()
df2['DormantN'] = df['DormantN'].groupby(df['sim']).mean()
df2['%Dormant'] = df2['DormantN']/df2['N']
df2['Prod'] = df['PRODI'].groupby(df['sim']).mean()
df2['ActiveN'] = df2['N'] - df2['DormantN']

#------------------------- DATA FILTERS -------------------

dat1 = pd.DataFrame(df2)

dat3 = df2[df2['R'] > 200]
dat3 = dat3[dat3['SpatialComplexity'].str.contains('-brownian-')]
dat3 = dat3[dat3['SpatialComplexity'].str.contains('-chemotaxis-')]
dat3 = dat3[dat3['ResourceComplexity'].str.contains('-lockandkey-')]

dat4 = df2[df2['R'] <= 200]
dat4 = dat4[dat4['SpatialComplexity'].str.contains('-brownian-')]
dat4 = dat4[dat4['SpatialComplexity'].str.contains('-chemotaxis-')]
dat4 = dat4[dat4['ResourceComplexity'].str.contains('-lockandkey-')]

#-------------------------------------------------------------------------------

#### plot figure ###############################################################
fs = 8 # fontsize
fig = plt.figure()
gd = 40
mct = 1
bns = 'log'

#### PLOT 1 #################################################################
fig.add_subplot(2, 2, 1)
xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
#ylab = '% Dormancy, '+'$log$'+r'$_{10}$'
ylab = '% Dormancy'


x1 = dat1['Encounters']
y1 = dat1['%Dormant']
plt.hexbin(x1, y1, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.winter)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
#plt.ylim(0.18, 0.7)
plt.ylim(0.1, 0.7)

plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.legend(bbox_to_anchor=(-0.04, 1.05, 2.48, .2), loc=10, ncol=3, mode="expand",prop={'size':fs})

#### PLOT 2 ################################
fig.add_subplot(2, 2, 2)

xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
#ylab = '% Dormancy, '+'$log$'+r'$_{10}$'
ylab = '% Dormancy'


x1 = dat3['Encounters']
y1 = dat3['%Dormant']
plt.hexbin(x1, y1, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.winter)

x2 = dat4['Encounters']
y2 = dat4['%Dormant']
plt.hexbin(x2, y2, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.autumn)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
#plt.ylim(0.18, 0.7)
plt.ylim(0.3, 0.7)

plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.legend(bbox_to_anchor=(-0.04, 1.05, 2.48, .2), loc=10, ncol=3, mode="expand",prop={'size':fs})

#### PLOT 3 #################################################################
fig.add_subplot(2, 2, 3)

xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
ylab = 'Total abundance, '+'$log$'+r'$_{10}$'

x1 = dat1['Encounters']
y1 = np.log10(dat1['N'])
plt.hexbin(x1, y1, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.winter)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.ylim(2.0, 3.4)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### PLOT 4 #################################################################
fig.add_subplot(2, 2, 4)

xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
ylab = 'Total abundance, '+'$log$'+r'$_{10}$'

x1 = dat3['Encounters']
y1 = np.log10(dat3['N'])
plt.hexbin(x1, y1, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.winter)

x2 = dat4['Encounters']
y2 = np.log10(dat4['N'])
plt.hexbin(x2, y2, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.autumn)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.ylim(2.0, 3.4)
plt.ylim(2.4, 3.4)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/EncounterFig_Heat.png', dpi=600, bbox_inches = "tight")
plt.close()
