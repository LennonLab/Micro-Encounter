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

dat = dat[dat['ResourceComplexityLevel'] != 3]
#dat = dat[dat['TrophicComplexityLevel'] == 4]
dat = dat[dat['SpatialComplexityLevel'] != 1]

#-------------------------------------------------------------------------------

#### plot figure ###############################################################
fs = 8 # fontsize
fig = plt.figure()
gd = 40
fs = 8

#### PLOT 1 #################################################################
fig.add_subplot(2, 2, 1)

ylab = 'Average encounters, '+'$log$'+r'$_{10}$'
xlab = 'Resource Inflow, '+'$log$'+r'$_{10}$'

y = np.log10(dat['MeanEncounter'])
x = np.log10(dat['ResInflow'])
plt.hexbin(x, y, mincnt=1, gridsize = gd, bins='log', cmap=plt.cm.jet)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
#plt.xlim(0.1, 300)
plt.tick_params(axis='both', which='major', labelsize=fs+2)

#### PLOT 2 #################################################################
fig.add_subplot(2, 2, 2)

ylab = 'Average encounters, '+'$log$'+r'$_{10}$'
xlab = 'Resource aggregation, '+'$log$'+r'$_{10}$'

y = np.log10(dat['MeanEncounter'])
x = np.log10(dat['MeanResAgg'])
plt.hexbin(x, y, mincnt=1, gridsize = gd, bins='log', cmap=plt.cm.jet)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.tick_params(axis='both', which='major', labelsize=fs+2)

#### PLOT 3 #################################################################
fig.add_subplot(2, 2, 3)

ylab = 'Average encounters, '+'$log$'+r'$_{10}$'
xlab = 'Total resources , '+'$log$'+r'$_{10}$'

y = np.log10(dat['MeanEncounter'])
x = np.log10(dat['MeanResourceParticles'])
plt.hexbin(x, y, mincnt=1, gridsize = gd, bins='log', cmap=plt.cm.jet)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.tick_params(axis='both', which='major', labelsize=fs+2)

#### PLOT 4 #################################################################
fig.add_subplot(2, 2, 4)

ylab = 'Average encounters, '+'$log$'+r'$_{10}$'
xlab = 'Incoming resource aggregation'

y = np.log10(dat['MeanEncounter'])
x = dat['IncomingResAgg']
plt.hexbin(x, y, mincnt=1, gridsize = gd, bins='log', cmap=plt.cm.jet)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.tick_params(axis='both', which='major', labelsize=fs+2)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/Encounter_vs_ResVars_No-SC1_RC3.png', dpi=600, bbox_inches = "tight")
#plt.show()
#plt.close()
