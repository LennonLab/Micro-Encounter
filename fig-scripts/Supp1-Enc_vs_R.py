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
#df = df[df['sim'] > 500]

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

df2['R'] = np.log10(df['R'].groupby(df['sim']).mean())
df2 = df2[np.isfinite(df2['R'])]

df2['ResInflow'] = df['ResInflow'].groupby(df['sim']).mean()
df2 = df2[np.isfinite(df2['ResInflow'])]

#------------------------- DATA FILTERS -------------------

dat0 = pd.DataFrame(df2)

#-------------------------------------------------------------------------------

#### plot figure ###############################################################
fs = 8 # fontsize
fig = plt.figure()
gd = 30
mct = 1
bns = 'log'

dat1 = df2[df2['SpatialComplexity'].str.contains('-wellmixed-')]
dat1 = dat1[dat1['SpatialComplexity'].str.contains('-none-')]
dat1 = dat1[dat1['ResourceComplexity'].str.contains('-lockandkey-')]

dat2 = df2[df2['SpatialComplexity'].str.contains('-wellmixed-')]
dat2 = dat2[dat2['SpatialComplexity'].str.contains('-chemotaxis-')]
dat2 = dat2[dat2['ResourceComplexity'].str.contains('-simple-')]

#### PLOT 1 #################################################################
fig.add_subplot(2, 2, 1)

#xlab = 'Resources, '+'$log$'+r'$_{10}$'
xlab = 'Resource inflow, '+'$log$'+r'$_{10}$'
ylab = 'Encounters, '+'$log$'+r'$_{10}$'

p1 = Rectangle((0, 0), 1, 1, fc='b', ec='b')
p2 = Rectangle((0, 0), 1, 1, fc='r', ec='r')
plt.legend([p1, p2], ['Recalcitrant + No dispersal', 'Labile + Chemotaxis'], bbox_to_anchor=(-0.03, 1.05, 2.46, .2), loc=10, ncol=2, mode="expand",prop={'size':fs+2})

#plt.xlim(0.8, 3.9)
plt.xlim(0.1, 0.9)
#xlim_list = [0.4, 0.4, 4.2, 4.2]
xlim_list = [0.0, 1.5, 0.0, 1.5]
plt.ylim(-1.71, 2.11)
ylim_list = [-1.75, 2.2, -1.75, 2.2]

#x0 = dat0['R'].tolist()
x0 = dat0['ResInflow'].tolist()
y0 = dat0['Encounters'].tolist()

#print min(x0), max(x0)
#print min(y0), max(y0)
#sys.exit()

x0.extend(xlim_list)
y0.extend(ylim_list)
plt.hexbin(x0, y0, mincnt=0, gridsize = gd, bins=bns, cmap=plt.cm.binary)

#x1 = dat1['R'].tolist()
x1 = dat1['ResInflow'].tolist()
y1 = dat1['Encounters'].tolist()

x1.extend(xlim_list)
y1.extend(ylim_list)
plt.hexbin(x1, y1, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.winter)

#x2 = dat2['R'].tolist()
x2 = dat2['ResInflow'].tolist()
y2 = dat2['Encounters'].tolist()

x2.extend(xlim_list)
y2.extend(ylim_list)
plt.hexbin(x2, y2, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.autumn)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)

plt.tick_params(axis='both', which='major', labelsize=fs)

#plt.text(1.6, -3.2, 'Well-mixed', fontsize=fs+10)
plt.text(0.3, -3.2, 'Well-mixed', fontsize=fs+10)

dat1 = df2[~df2['SpatialComplexity'].str.contains('-wellmixed-')]
dat1 = dat1[dat1['SpatialComplexity'].str.contains('-none-')]
dat1 = dat1[dat1['ResourceComplexity'].str.contains('-lockandkey-')]

dat2 = df2[~df2['SpatialComplexity'].str.contains('-wellmixed-')]
dat2 = dat2[dat2['SpatialComplexity'].str.contains('-chemotaxis-')]
dat2 = dat2[dat2['ResourceComplexity'].str.contains('-simple-')]


#### PLOT 2 #################################################################
fig.add_subplot(2, 2, 2)
#xlab = 'Resources, '+'$log$'+r'$_{10}$'
xlab = 'Resource inflow, '+'$log$'+r'$_{10}$'
ylab = 'Encounters, '+'$log$'+r'$_{10}$'

p1 = Rectangle((0, 0), 1, 1, fc='b', ec='b')
p2 = Rectangle((0, 0), 1, 1, fc='r', ec='r')

#plt.xlim(0.8, 3.9)
plt.xlim(0.1, 0.9)
#xlim_list = [0.4, 0.4, 4.2, 4.2]
xlim_list = [0.0, 1.5, 0.0, 1.5]
plt.ylim(-1.71, 2.11)
ylim_list = [-1.75, 2.2, -1.75, 2.2]

#x0 = dat0['R'].tolist()
x0 = dat0['ResInflow'].tolist()
y0 = dat0['Encounters'].tolist()

#print min(x0), max(x0)
#print min(y0), max(y0)
#sys.exit()

x0.extend(xlim_list)
y0.extend(ylim_list)
plt.hexbin(x0, y0, mincnt=0, gridsize = gd, bins=bns, cmap=plt.cm.binary)

#x1 = dat1['R'].tolist()
x1 = dat1['ResInflow'].tolist()
y1 = dat1['Encounters'].tolist()

x1.extend(xlim_list)
y1.extend(ylim_list)
plt.hexbin(x1, y1, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.winter)

#x2 = dat2['R'].tolist()
x2 = dat2['ResInflow'].tolist()
y2 = dat2['Encounters'].tolist()

x2.extend(xlim_list)
y2.extend(ylim_list)
plt.hexbin(x2, y2, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.autumn)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)

plt.tick_params(axis='both', which='major', labelsize=fs)

#plt.text(1.6, -3.2, 'Structured', fontsize=fs+10)
plt.text(0.3, -3.2, 'Structured', fontsize=fs+10)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/Supp1-Env_vs_ResInflow.png', dpi=600, bbox_inches = "tight")
plt.close()
