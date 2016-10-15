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
df2 = df2[np.isfinite(df2['D'])]

df2['N'] = df['N'].groupby(df['sim']).mean()
df2 = df2[np.isfinite(df2['N'])]

df2['DormantN'] = df['DormantN'].groupby(df['sim']).mean()
df2 = df2[np.isfinite(df2['DormantN'])]

df2['%Dormant'] = df2['DormantN']/df2['N']
df2 = df2[np.isfinite(df2['%Dormant'])]

df2['Prod'] = df['PRODI'].groupby(df['sim']).mean()
df2 = df2[np.isfinite(df2['Prod'])]

df2['ActiveN'] = df2['N'] - df2['DormantN']
df2 = df2[np.isfinite(df2['ActiveN'])]

df2['Q'] = df['MeanCellQuota'].groupby(df['sim']).mean()
df2 = df2[np.isfinite(df2['Q'])]

df2['R'] = df['R'].groupby(df['sim']).mean()
df2 = df2[np.isfinite(df2['R'])]

#------------------------- DATA FILTERS -------------------

dat0 = pd.DataFrame(df2)
df2 = df2[~df2['SpatialComplexity'].str.contains('-wellmixed-')]

dat1 = df2[df2['SpatialComplexity'].str.contains('-none-')]
dat1 = dat1[dat1['ResourceComplexity'].str.contains('-lockandkey-')]

dat2 = df2[df2['SpatialComplexity'].str.contains('-chemotaxis-')]
#dat2 = dat2[dat2['SpatialComplexity'].str.contains('-brownian-')]
dat2 = dat2[dat2['ResourceComplexity'].str.contains('-simple-')]
#dat2 = dat2[dat2['TrophicComplexity'].str.contains('-scavenging-')]

#-------------------------------------------------------------------------------

#### plot figure ###############################################################
fs = 8 # fontsize
fig = plt.figure()
gd = 30
mct = 1
bns = 'log'

#### PLOT 1 #################################################################
fig.add_subplot(2, 2, 1)
xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
#ylab = '% Dormancy, '+'$log$'+r'$_{10}$'
ylab = '% Dormant'

p1 = Rectangle((0, 0), 1, 1, fc='b', ec='b')
p2 = Rectangle((0, 0), 1, 1, fc='r', ec='r')
plt.legend([p1, p2], ['Recalcitrant + No dispersal', 'Labile + Chemotaxis'], bbox_to_anchor=(-0.04, 1.05, 2.48, .2), loc=10, ncol=2, mode="expand",prop={'size':fs+5})

plt.xlim(-1.6, 2.2)
xlim_list = [-1.8, -1.8, 2.4, 2.4]
plt.ylim(0.0, 0.6)
ylim_list = [-0.1, 0.6, -0.1, 0.6]

x0 = dat0['Encounters'].tolist()
y0 = dat0['%Dormant'].tolist()
x0.extend(xlim_list)
y0.extend(ylim_list)
plt.hexbin(x0, y0, mincnt=0, gridsize = gd, bins=bns, cmap=plt.cm.binary)

x1 = dat1['Encounters'].tolist()
y1 = dat1['%Dormant'].tolist()
x1.extend(xlim_list)
y1.extend(ylim_list)
plt.hexbin(x1, y1, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.winter)

x2 = dat2['Encounters'].tolist()
y2 = dat2['%Dormant'].tolist()
x2.extend(xlim_list)
y2.extend(ylim_list)
plt.hexbin(x2, y2, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.autumn)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)

plt.tick_params(axis='both', which='major', labelsize=fs)

plt.text(-0.5, -0.23, 'Microscale', fontsize=fs+12)

#### PLOT 2 ################################
fig.add_subplot(2, 2, 2)

ylab = '% Dormant'
xlab = 'Resources, '+'$log$'+r'$_{10}$'

plt.xlim(0.8, 4.0)
xlim_list = [0.8, 0.8, 4.2, 4.2]
plt.ylim(0.0, 0.61)
ylim_list = [-0.1, 0.64, -0.1, 0.64]

y0 = dat0['%Dormant'].tolist()
x0 = np.log10(dat0['R']).tolist()

#print min(x0), max(x0)
#print min(y0), max(y0)
#sys.exit()

x0.extend(xlim_list)
y0.extend(ylim_list)
plt.hexbin(x0, y0, mincnt=0, gridsize = gd, bins=bns, cmap=plt.cm.binary)

y1 = dat1['%Dormant'].tolist()
x1 = np.log10(dat1['R']).tolist()
x1.extend(xlim_list)
y1.extend(ylim_list)
plt.hexbin(x1, y1, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.winter)

y2 = dat2['%Dormant'].tolist()
x2 = np.log10(dat2['R']).tolist()
x2.extend(xlim_list)
y2.extend(ylim_list)
plt.hexbin(x2, y2, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.autumn)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.tick_params(axis='both', which='major', labelsize=fs)

plt.text(1.75, -0.23, 'Macroscale', fontsize=fs+12)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/MicroMacro_Heat.png', dpi=600, bbox_inches = "tight")
plt.close()
