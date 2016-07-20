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
#-------------------------DATA FILTERS------------------------------------------

#mdat = mdat[mdat['ResInflow'] > 20]
#mdat = mdat[mdat['MaxMetMaint'] < 0.002]
#mdat = mdat[mdat['IncomingResAgg'] < 0.4]
#mdat = mdat[mdat['MeanTotalAbundance'] > 200]

#-------------------------END DATA FILTERS--------------------------------------

color1 = 'm'
color2 = 'steelblue'
color3 = 'goldenrod'

#### plot figure ###############################################################
fs = 6 # fontsize
fig = plt.figure()

ax = fig.add_subplot(3, 3, 1)



dat = mdat[mdat['ResourceComplexityLevel'] == 3]
dat = dat[dat['TrophicComplexityLevel'] >= 3]

dat['Encounter'] = np.log10(dat['MeanEncounter'])
dat = dat[np.isfinite(dat['Encounter'])]

dat1 = dat[dat['SpatialComplexityLevel'] == 1]
encounters1 = dat1['Encounter']
dat2 = dat[dat['SpatialComplexityLevel'] == 2]
encounters2 = dat2['Encounter']
dat3 = dat[dat['SpatialComplexityLevel'] == 3]
encounters3 = dat3['Encounter']

data_to_plot = [encounters1, encounters2, encounters3]
ax.boxplot(data_to_plot, showmeans=True, showfliers=False)

ax.set_xticklabels(['White noise', 'Aggregated Resources +\nRandom walk', 'Aggregated Resources +\nChemotaxis'], rotation=45)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()


ax = fig.add_subplot(3, 3, 2)

dat = mdat[mdat['SpatialComplexityLevel'] == 3]
dat = dat[dat['TrophicComplexityLevel'] >= 3]

dat['Encounter'] = np.log10(dat['MeanEncounter'])
dat = dat[np.isfinite(dat['Encounter'])]

dat1 = dat[dat['ResourceComplexityLevel'] == 1]
encounters1 = dat1['Encounter']
dat2 = dat[dat['ResourceComplexityLevel'] == 2]
encounters2 = dat2['Encounter']
dat3 = dat[dat['ResourceComplexityLevel'] == 3]
encounters3 = dat3['Encounter']

data_to_plot = [encounters1, encounters2, encounters3]
ax.boxplot(data_to_plot, showmeans=True, showfliers=False)
ax.set_xticklabels(['Monoculture', 'Polyculture', 'Lock & Key'], rotation=45)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()



ax = fig.add_subplot(3, 3, 3)

dat = mdat[mdat['SpatialComplexityLevel'] <= 3]
dat = dat[dat['ResourceComplexityLevel'] <= 2] 

dat['Encounter'] = np.log10(dat['MeanEncounter'])
dat = dat[np.isfinite(dat['Encounter'])]

dat1 = dat[dat['TrophicComplexityLevel'] == 1]
encounters1 = dat1['Encounter']
dat2 = dat[dat['TrophicComplexityLevel'] == 2]
encounters2 = dat2['Encounter']
dat3 = dat[dat['TrophicComplexityLevel'] == 3]
encounters3 = dat3['Encounter']
dat4 = dat[dat['TrophicComplexityLevel'] == 4]
encounters4 = dat4['Encounter']

data_to_plot = [encounters1, encounters2, encounters3, encounters4]
ax.boxplot(data_to_plot, showmeans=True, showfliers=False)
ax.set_xticklabels(['Consumer-Resource', 'Crossfeeding', 'Scavenging', 'Crossfeeding +\nScavenging'], rotation=45)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()


#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/BoxPlot.png', dpi=600, bbox_inches = "tight")
#plt.show()
#plt.close()
