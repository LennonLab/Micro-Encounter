from __future__ import division
import  matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd
import numpy as np
import os
import sys

mydir = os.path.expanduser('~/GitHub/Micro-Encounter')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")
#df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')
df = pd.read_csv(mydir + '/results/simulated_data/2016_09_18_SimData.csv')
#df = df.convert_objects(convert_numeric=True).dropna()
df = df[df['sim'] > 500]

#-------------------------DATA TRANSFORMATIONS -----------------------
df2 = pd.DataFrame({'Encounters' : np.log10(df['Encounters'].groupby(df['sim']).mean())})
df2 = df2[np.isfinite(df2['Encounters'])]

df2['SpatialComplexity'] = df['SpatialComplexity'].groupby(df['sim']).unique()
df2['SpatialComplexity'] = [df2['SpatialComplexity'][i][0] for i in df2['SpatialComplexity'].keys()]

df2['TrophicComplexity'] = df['TrophicComplexity'].groupby(df['sim']).unique()
df2['TrophicComplexity'] = [df2['TrophicComplexity'][i][0] for i in df2['TrophicComplexity'].keys()]

df2['ResourceComplexity'] = df['ResourceComplexity'].groupby(df['sim']).unique()
df2['ResourceComplexity'] = [df2['ResourceComplexity'][i][0] for i in df2['ResourceComplexity'].keys()]

df2['D'] = df['numDead'].groupby(df['sim']).mean()
df2['N'] = df['N'].groupby(df['sim']).mean()
df2['DormantN'] = df['DormantN'].groupby(df['sim']).mean()
df2['%Dormant'] = df2['DormantN']/df2['N']
df2['Prod'] = df['PRODI'].groupby(df['sim']).mean()
df2['ActiveN'] = df2['N'] - df2['DormantN']
df2['R'] = df['R'].groupby(df['sim']).mean()
df2['ResInflow'] = df['ResInflow'].groupby(df['sim']).mean()

# TRAITS
df2['Q'] = df['MeanCellQuota'].groupby(df['sim']).mean()
df2['MaxGrowth'] = df['MaxGrowth'].groupby(df['sim']).mean()
df2['MaxMaint'] = df['MaxMaint'].groupby(df['sim']).mean()
df2['MaxDispersal'] = df['MaxDispersal'].groupby(df['sim']).mean()
df2['MaxRPF'] = df['MaxRPF'].groupby(df['sim']).mean()
df2['MaxMainFactor'] = df['MaxMainFactor'].groupby(df['sim']).mean()
df2['SpecDisp'] = df['SpeciesDisp'].groupby(df['sim']).mean()
df2['SpecMaint'] = df['SpeciesMaint'].groupby(df['sim']).mean()
df2['SpecGrowth'] = df['SpecificGrowth'].groupby(df['sim']).mean()
df2['PerCapitaGrowth'] = df['PerCapitaGrowth'].groupby(df['sim']).mean()
df2['PerCapitaMaint'] = df['PerCapitaMaint'].groupby(df['sim']).mean()
df2['PerCapitaDisp'] = df['PerCapitaDisp'].groupby(df['sim']).mean()
df2['numDead'] = df['numDead'].groupby(df['sim']).mean()

#------------------------- DATA FILTERS -------------------
dat0 = pd.DataFrame(df2)

df2 = df2[df2['SpatialComplexity'].str.contains('-wellmixed-')]

dat1 = df2[df2['SpatialComplexity'].str.contains('-none-')]
#dat1 = dat1[dat1['ResourceComplexity'].str.contains('-lockandkey-')]

dat2 = df2[df2['SpatialComplexity'].str.contains('-chemotaxis-')]
#dat2 = dat2[dat2['SpatialComplexity'].str.contains('-brownian-')]
#dat2 = dat2[dat2['ResourceComplexity'].str.contains('-simple-')]
#dat2 = dat2[dat2['TrophicComplexity'].str.contains('-scavenging-')]

#-------------------------------------------------------------------------------

#### plot figure ###############################################################
fs = 4 # fontsize
fig = plt.figure()
gd = 30
mct = 1
bns = 'log'

#### PLOT 1 #################################################################
fig.add_subplot(3, 3, 1)
xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
#ylab = '% Dormancy, '+'$log$'+r'$_{10}$'
ylab = 'Dispersal'

p1 = Rectangle((0, 0), 1, 1, fc='b', ec='b')
p2 = Rectangle((0, 0), 1, 1, fc='r', ec='r')
plt.legend([p1, p2], ['Recalcitrant + No dispersal', 'Labile + Chemotaxis'], bbox_to_anchor=(-0.04, 1.05, 3.88, .2), loc=10, ncol=2, mode="expand",prop={'size':fs+5})

plt.xlim(-1.75, 2.15)
xlim_list = [-1.95, -1.95, 2.35, 2.35]
plt.ylim(-0.11, 0.9)
ylim_list = [-0.15, 0.99, -0.15, 0.99]

x0 = dat0['Encounters'].tolist()
y0 = dat0['PerCapitaDisp'].tolist()
#y0 = dat0['SpecDisp'].tolist()

x0.extend(xlim_list)
y0.extend(ylim_list)
plt.hexbin(x0, y0, mincnt=0, gridsize = gd, bins=bns, cmap=plt.cm.binary)

x1 = dat1['Encounters'].tolist()
y1 = dat1['PerCapitaDisp'].tolist()
#y1 = dat1['SpecDisp'].tolist()

x1.extend(xlim_list)
y1.extend(ylim_list)
plt.hexbin(x1, y1, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.winter)

x2 = dat2['Encounters'].tolist()
y2 = dat2['PerCapitaDisp'].tolist()
#y2 = dat2['SpecDisp'].tolist()

x2.extend(xlim_list)
y2.extend(ylim_list)
plt.hexbin(x2, y2, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.autumn)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)

plt.tick_params(axis='both', which='major', labelsize=fs)
#plt.legend(bbox_to_anchor=(-0.04, 1.05, 2.48, .2), loc=10, ncol=3, mode="expand",prop={'size':fs})

plt.text(-3.3, 0.7, 'Per capita', rotation = 90, fontsize=fs+8)

#### PLOT 2 ################################
fig.add_subplot(3, 3, 2)

xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
ylab = 'Maintenance'

plt.xlim(-2.0, 2.15)
xlim_list = [-2.95, -2.95, 2.35, 2.35]
plt.ylim(0.7, 6.8)
ylim_list = [0.5, 7.0, 0.5, 7.0]

x0 = dat0['Encounters'].tolist()
y0 = dat0['PerCapitaMaint'].tolist()
#y0 = dat0['SpecMaint'].tolist()

x0.extend(xlim_list)
y0.extend(ylim_list)
plt.hexbin(x0, y0, mincnt=0, gridsize = gd, bins=bns, cmap=plt.cm.binary)

x1 = dat1['Encounters'].tolist()
y1 = dat1['PerCapitaMaint'].tolist()
#y1 = dat1['SpecMaint'].tolist()

x1.extend(xlim_list)
y1.extend(ylim_list)
plt.hexbin(x1, y1, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.winter)

x2 = dat2['Encounters'].tolist()
y2 = dat2['PerCapitaMaint'].tolist()
#y2 = dat2['SpecMaint'].tolist()

x2.extend(xlim_list)
y2.extend(ylim_list)
plt.hexbin(x2, y2, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.autumn)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.tick_params(axis='both', which='major', labelsize=fs)


#### PLOT 3 #################################################################
fig.add_subplot(3, 3, 3)

xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
ylab = 'Growth'

plt.xlim(-2.0, 2.15)
xlim_list = [-2.95, -2.95, 2.35, 2.35]
plt.ylim(0.06, 0.9)
ylim_list = [0.04, 1.0, 0.04, 1.0]

x0 = dat0['Encounters'].tolist()
y0 = dat0['PerCapitaGrowth'].tolist()
#y0 = dat0['SpecGrowth'].tolist()

x0.extend(xlim_list)
y0.extend(ylim_list)

plt.hexbin(x0, y0, mincnt=0, gridsize = gd, bins=bns, cmap=plt.cm.binary)

x1 = dat1['Encounters'].tolist()
y1 = dat1['PerCapitaGrowth'].tolist()
#y1 = dat1['SpecGrowth'].tolist()

x1.extend(xlim_list)
y1.extend(ylim_list)
plt.hexbin(x1, y1, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.winter)

x2 = dat2['Encounters'].tolist()
y2 = dat2['PerCapitaGrowth'].tolist()
#y2 = dat2['SpecGrowth'].tolist()

x2.extend(xlim_list)
y2.extend(ylim_list)
plt.hexbin(x2, y2, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.autumn)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### SECOND ROW #####

#### PLOT 1 #################################################################
fig.add_subplot(3, 3, 4)
xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
#ylab = '% Dormancy, '+'$log$'+r'$_{10}$'
ylab = 'Dispersal'

p1 = Rectangle((0, 0), 1, 1, fc='b', ec='b')
p2 = Rectangle((0, 0), 1, 1, fc='r', ec='r')

plt.xlim(-1.75, 2.15)
xlim_list = [-1.95, -1.95, 2.35, 2.35]
plt.ylim(-0.11, 0.9)
ylim_list = [-0.15, 0.99, -0.15, 0.99]

x0 = dat0['Encounters'].tolist()
#y0 = dat0['PerCapitaDisp'].tolist()
y0 = dat0['SpecDisp'].tolist()

x0.extend(xlim_list)
y0.extend(ylim_list)
plt.hexbin(x0, y0, mincnt=0, gridsize = gd, bins=bns, cmap=plt.cm.binary)

x1 = dat1['Encounters'].tolist()
#y1 = dat1['PerCapitaDisp'].tolist()
y1 = dat1['SpecDisp'].tolist()

x1.extend(xlim_list)
y1.extend(ylim_list)
plt.hexbin(x1, y1, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.winter)

x2 = dat2['Encounters'].tolist()
#y2 = dat2['PerCapitaDisp'].tolist()
y2 = dat2['SpecDisp'].tolist()

x2.extend(xlim_list)
y2.extend(ylim_list)
plt.hexbin(x2, y2, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.autumn)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)

plt.tick_params(axis='both', which='major', labelsize=fs)

plt.text(-3.3, 0.8, 'Species specific', rotation = 90, fontsize=fs+8)

#### PLOT 2 ################################
fig.add_subplot(3, 3, 5)

xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
ylab = 'Maintenance'

plt.xlim(-2.0, 2.15)
xlim_list = [-2.95, -2.95, 2.35, 2.35]
plt.ylim(0.7, 6.8)
ylim_list = [0.5, 7.0, 0.5, 7.0]

x0 = dat0['Encounters'].tolist()
#y0 = dat0['PerCapitaMaint'].tolist()
y0 = dat0['SpecMaint'].tolist()

x0.extend(xlim_list)
y0.extend(ylim_list)
plt.hexbin(x0, y0, mincnt=0, gridsize = gd, bins=bns, cmap=plt.cm.binary)

x1 = dat1['Encounters'].tolist()
#y1 = dat1['PerCapitaMaint'].tolist()
y1 = dat1['SpecMaint'].tolist()

x1.extend(xlim_list)
y1.extend(ylim_list)
plt.hexbin(x1, y1, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.winter)

x2 = dat2['Encounters'].tolist()
#y2 = dat2['PerCapitaMaint'].tolist()
y2 = dat2['SpecMaint'].tolist()

x2.extend(xlim_list)
y2.extend(ylim_list)
plt.hexbin(x2, y2, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.autumn)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.tick_params(axis='both', which='major', labelsize=fs)


#### PLOT 3 #################################################################
fig.add_subplot(3, 3, 6)

xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
ylab = 'Growth'

plt.xlim(-2.0, 2.15)
xlim_list = [-2.95, -2.95, 2.35, 2.35]
plt.ylim(0.06, 0.9)
ylim_list = [0.04, 1.0, 0.04, 1.0]

x0 = dat0['Encounters'].tolist()
#y0 = dat0['PerCapitaGrowth'].tolist()
y0 = dat0['SpecGrowth'].tolist()

x0.extend(xlim_list)
y0.extend(ylim_list)

plt.hexbin(x0, y0, mincnt=0, gridsize = gd, bins=bns, cmap=plt.cm.binary)

x1 = dat1['Encounters'].tolist()
#y1 = dat1['PerCapitaGrowth'].tolist()
y1 = dat1['SpecGrowth'].tolist()

x1.extend(xlim_list)
y1.extend(ylim_list)
plt.hexbin(x1, y1, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.winter)

x2 = dat2['Encounters'].tolist()
#y2 = dat2['PerCapitaGrowth'].tolist()
y2 = dat2['SpecGrowth'].tolist()

x2.extend(xlim_list)
y2.extend(ylim_list)
plt.hexbin(x2, y2, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.autumn)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/Traits_Heat.png', dpi=600, bbox_inches = "tight")
plt.close()
