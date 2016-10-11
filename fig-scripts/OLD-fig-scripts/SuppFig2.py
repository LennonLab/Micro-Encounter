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
mdat = mdat[mdat['RowID'] > 3000]
#print mdat.shape



#### plot figure ###############################################################
dat = mdat[mdat['ResourceComplexityLevel'] == 3] # Chemical complexity
dat = dat[dat['SpatialComplexityLevel'] == 3] # Chemotaxis


fs = 6 # fontsize
fig = plt.figure()
ax = fig.add_subplot(2, 2, 1)

y = dat['MeanPerCapitaActiveDispersal']
x = np.log10(dat['ResInflow'])

plt.scatter(x, y, color = '0.2', alpha = 0.2, label='chemotaxis')
plt.ylabel('Per Capita Dispersal')
plt.xlabel('Resource Inflow')
#plt.xlim(-2, 1)
plt.tick_params(axis='both', which='major', labelsize=fs+4)



ax = fig.add_subplot(2, 2, 2)

y = dat['MeanPerCapitaActiveDispersal']
x = dat['MeanIndProduction']

plt.scatter(x, y, color = '0.2', alpha = 0.2, label='chemotaxis')
plt.ylabel('Per Capita Dispersal')
plt.xlabel('Productivity')
plt.xlim(-0.2, 0.6)
plt.tick_params(axis='both', which='major', labelsize=fs+4)



ax = fig.add_subplot(2, 2, 3)

y = dat['MeanSpecDisp']
x = np.log10(dat['ResInflow'])

plt.scatter(x, y, color = '0.2', alpha = 0.2, label='chemotaxis')
plt.ylabel('Species Specific Dispersal')
plt.xlabel('Resource Inflow')
plt.tick_params(axis='both', which='major', labelsize=fs+4)




ax = fig.add_subplot(2, 2, 4)

y = dat['MeanSpecDisp']
x = dat['MeanIndProduction']

plt.scatter(x, y, color = '0.2', alpha = 0.2, label='chemotaxis')
plt.ylabel('Species Specific Dispersal')
plt.xlabel('Productivity')
plt.xlim(-0.2, 0.6)
plt.tick_params(axis='both', which='major', labelsize=fs+4)


#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.45, hspace=0.4)
plt.savefig(mydir + '/results/figures/SupplementalFigure1.png', dpi=600, bbox_inches = "tight")
#plt.show()
#plt.close()
