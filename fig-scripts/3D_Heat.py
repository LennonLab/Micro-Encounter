from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

from statsmodels.nonparametric.smoothers_lowess import lowess


import numpy as np

import pandas as pd
import os
import sys

import scipy as sc
from scipy import stats

mydir = os.path.expanduser('~/GitHub/Micro-Encounter')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")
#dat = pd.read_csv(mydir + '/results/simulated_data/2016_07_19_SimData.csv')
dat = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')
dat = dat.convert_objects(convert_numeric=True).dropna()


#-------------------------DATA TRANSFORMS---------------------------------------

dat['DormFreq'] = np.log10(dat['MeanDormFreq'])
dat = dat[np.isfinite(dat['DormFreq'])]

dat['Encounters'] = np.log10(dat['VarEncounter'])
dat = dat[np.isfinite(dat['Encounters'])]
dat = dat[dat['Encounters'] < 3]

dat['Res_Inflow'] = np.log10(dat['ResInflow'])
dat = dat[np.isfinite(dat['Res_Inflow'])]

dat['AvgResourceParticles'] = np.log10(dat['MeanResourceParticles'])
dat = dat[np.isfinite(dat['AvgResourceParticles'])]

dat['Production'] = np.log10(dat['MeanIndProduction'])
dat = dat[np.isfinite(dat['Production'])]

dat['TotalAbundance'] = np.log10(dat['MeanTotalAbundance'])
dat = dat[np.isfinite(dat['TotalAbundance'])]

dat['ActiveAbundance'] = np.log10(dat['MeanTotalAbundance'] * (1 - dat['MeanDormFreq']))
dat = dat[np.isfinite(dat['ActiveAbundance'])]

#-------------------------END DATA TRANSFORMS-----------------------------------

#-------------------------DATA FILTERS------------------------------------------

#dat = dat[dat['ResourceComplexityLevel'] == 3]
#dat = dat[dat['TrophicComplexityLevel'] == 4]
#dat = dat[dat['SpatialComplexityLevel'] == 2]
print 'size of dat:', dat.shape

#-------------------------END DATA FILTERS--------------------------------------


#### figure ###############################################################
fig = plt.figure()
ax = fig.gca(projection='3d')

x = dat['Res_Inflow']
y = dat['Encounters']
z = dat['DormFreq']

#zfilt = lowess(x, z, is_sorted=False, frac=0.1, it=0)
#z = zfilt[:,0]

ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.0, antialiased=True)
plt.show()

#fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_xlabel('ResInflow')
ax.set_ylabel('Encounter')
ax.set_zlabel('%Dormancy')

plt.show()
#plt.close()
