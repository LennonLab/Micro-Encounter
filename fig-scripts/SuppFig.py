from __future__ import division
import  matplotlib.pyplot as plt
import random
from random import shuffle
import pandas as pd
import numpy as np
import os
import sys

import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import summary_table

mydir = os.path.expanduser('~/GitHub/Micro-Encounter')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")
mdat = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')
mdat = mdat.convert_objects(convert_numeric=True).dropna()

print mdat.shape


#### plot figure ###############################################################
fs = 6 # fontsize
fig = plt.figure()

ax = fig.add_subplot(2, 2, 1)

dat = mdat[mdat['ResourceComplexityLevel'] == 3]
dat = dat[dat['MeanResourceParticles'] < 150]
dat1 = dat[dat['SpatialComplexityLevel'] == 3]
dat1['Encounter'] = np.log10(dat1['MeanEncounter'])
dat1 = dat1[np.isfinite(dat1['Encounter'])]
dat1 = dat1[dat1['Encounter'] < 0.5]
y1 = dat1['Encounter']
x1 = np.log10(dat1['MeanResourceParticles'])

dat2 = dat[dat['SpatialComplexityLevel'] == 1]
dat2['Encounter'] = np.log10(dat2['MeanEncounter'])
dat2 = dat2[np.isfinite(dat2['Encounter'])]
dat2 = dat2[dat2['Encounter'] < 0.5]
y2 = dat2['Encounter']
x2 = np.log10(dat2['MeanResourceParticles'])

plt.scatter(x1, y1, color = 'm', alpha = 0.6, label='chemotaxis')
plt.scatter(x2, y2, color = '0.2', alpha = 0.6, label='white noise')
#plt.ylim(-1.5, 0.5)
plt.ylabel('Encounters')
plt.xlabel('Resources')
plt.tick_params(axis='both', which='major', labelsize=fs+4)

plt.legend(bbox_to_anchor=(-0.04, 1.05, 2.48, .2), loc=10, ncol=2, mode="expand",prop={'size':fs+6})


ax = fig.add_subplot(2, 2, 2)

dat = mdat[mdat['ResourceComplexityLevel'] == 3]
#dat = dat[dat['MeanResourceParticles'] < 70]
dat = dat[dat['ResInflow'] <= 2]
dat1 = dat[dat['SpatialComplexityLevel'] == 3]
dat2 = dat[dat['SpatialComplexityLevel'] == 1]

dat1 = dat[dat['SpatialComplexityLevel'] == 3]
dat1['Encounter'] = np.log10(dat1['MeanEncounter'])
dat1 = dat1[np.isfinite(dat1['Encounter'])]
dat1 = dat1[dat1['Encounter'] < 0.5]
y1 = dat1['Encounter']
x1 = np.log10(dat1['MeanResourceParticles'])

dat2 = dat[dat['SpatialComplexityLevel'] == 1]
dat2['Encounter'] = np.log10(dat2['MeanEncounter'])
dat2 = dat2[np.isfinite(dat2['Encounter'])]
dat2 = dat2[dat2['Encounter'] < 0.5]
y2 = dat2['Encounter']
x2 = np.log10(dat2['MeanResourceParticles'])

plt.scatter(x1, y1, color = 'm', alpha = 0.6, label='chemotaxis')
plt.scatter(x2, y2, color = '0.2', alpha = 0.6, label='white noise')
#plt.ylim(-1.5, 0.5)
plt.ylabel('Encounters')
plt.xlabel('Resources')

#ax.set_xticklabels(['White noise', 'Aggregated Resources +\nRandom walk', 'Aggregated Resources +\nChemotaxis'], rotation=45)
plt.tick_params(axis='both', which='major', labelsize=fs+4)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()


#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
#plt.savefig(mydir + '/results/figures/BoxPlot.png', dpi=600, bbox_inches = "tight")
plt.show()
#plt.close()
