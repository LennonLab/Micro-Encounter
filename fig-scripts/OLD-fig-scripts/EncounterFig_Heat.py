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

dat = dat[dat['ResourceComplexityLevel'] != 3]
#dat = dat[dat['TrophicComplexityLevel'] != 3]
#dat = dat[dat['SpatialComplexityLevel'] == 3]

#dat = dat[dat['height'] <= 6]
#dat = dat[dat['width'] <= 6]

dat['DormFreq'] = np.log10(dat['MeanDormFreq'])
#dat['DormFreq'] = dat['MeanDormFreq']
dat = dat[np.isfinite(dat['DormFreq'])]
dat = dat[dat['DormFreq'] > -2]

dat['Encounters'] = np.log10(dat['MeanEncounter'])
dat = dat[np.isfinite(dat['Encounters'])]
dat = dat[dat['Encounters'] < 3]

dat['Production'] = np.log10(dat['MeanIndProduction'])
dat = dat[np.isfinite(dat['Production'])]

dat['TotalAbundance'] = np.log10(dat['MeanTotalAbundance'])
dat = dat[np.isfinite(dat['TotalAbundance'])]

dat['ActiveAbundance'] = np.log10(dat['MeanTotalAbundance'] * (1 - dat['MeanDormFreq']))
dat = dat[np.isfinite(dat['ActiveAbundance'])]

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
ylab = '% Dormancy, '+'$log$'+r'$_{10}$'
#ylab = '% Dormancy'
width = 1

x = dat['Encounters'].tolist()
y = dat['DormFreq'].tolist()
plt.hexbin(x, y, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.jet)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
#plt.ylim(-1.0, 0.1)
plt.tick_params(axis='both', which='major', labelsize=fs)
plt.legend(bbox_to_anchor=(-0.04, 1.05, 2.48, .2), loc=10, ncol=3, mode="expand",prop={'size':fs})

#### PLOT 2 ################################
fig.add_subplot(2, 2, 2)

xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
ylab = 'Productivity, '+'$log$'+r'$_{10}$'

x = dat['Encounters'].tolist()
y = dat['Production'].tolist()
plt.hexbin(x, y, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.jet)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
#plt.ylim(-2.0, 2.0)
#plt.xlim(0.1, 300)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### PLOT 3 #################################################################
fig.add_subplot(2, 2, 3)

xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
ylab = 'Total abundance, '+'$log$'+r'$_{10}$'

x = dat['Encounters'].tolist()
y = dat['TotalAbundance'].tolist()
plt.hexbin(x, y, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.jet)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
#plt.xlim(0.15, 300)
#plt.ylim(0.5, 3.1)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### PLOT 4 #################################################################
fig.add_subplot(2, 2, 4)

xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
ylab = 'Active abundance, '+'$log$'+r'$_{10}$'

x = dat['Encounters'].tolist()
y = dat['ActiveAbundance'].tolist()
plt.hexbin(x, y, mincnt=mct, gridsize = gd, bins=bns, cmap=plt.cm.jet)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
#plt.xlim(0.15, 1000)
#plt.ylim(-0.5, 3.1)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/EncounterFig_Heat_Spatial_RC2-SC3.png', dpi=600, bbox_inches = "tight")

#plt.show()
#plt.close()
