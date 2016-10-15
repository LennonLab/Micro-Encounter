from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
import numpy as np
import pandas as pd
import scipy as sc
from scipy import stats
import sys
import os


mydir = os.path.expanduser('~/GitHub/Micro-Encounter')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")
#dat = pd.read_csv(mydir + '/results/simulated_data/2016_07_19_SimData.csv')
#df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')
#df = pd.read_csv(mydir + '/results/simulated_data/SimData-TSA-Extremes.csv')
df = pd.read_csv(mydir + '/results/simulated_data/2016_09_18_SimData.csv')

#-------------------------DATA TRANSFORMS---------------------------------------

df2 = pd.DataFrame({'Encounters' : np.log10(df['Encounters'].groupby(df['sim']).mean())})
df2 = df2[np.isfinite(df2['Encounters'])]

#df2['R'] = np.log10(df['R'].groupby(df['sim']).mean())
#df2 = df2[np.isfinite(df2['R'])]

df2['R'] = np.log10(df['R'].groupby(df['sim']).mean())
df2 = df2[np.isfinite(df2['R'])]

df2['ResIn'] = np.log10(df['ResInflow'].groupby(df['sim']).mean())
#print 10**df2['ResIn']
#sys.exit()

#df2['P'] = np.log10(df['PRODI'].groupby(df['sim']).mean())
#df2 = df2[np.isfinite(df2['P'])]


df2['DormantN'] = df['DormantN'].groupby(df['sim']).mean()
df2 = df2[np.isfinite(df2['DormantN'])]

df2['N'] = df['N'].groupby(df['sim']).mean()
df2 = df2[np.isfinite(df2['N'])]

df2['DormFreq'] = np.log10(df2['DormantN']/df2['N'])
df2 = df2[np.isfinite(df2['DormFreq'])]

df2['N'] = np.log10(df2['N'])
df2 = df2[np.isfinite(df2['N'])]

#-------------------------END DATA TRANSFORMS-----------------------------------

#-------------------------DATA FILTERS------------------------------------------
print 'size of dat:', df2.shape

#print df2
#sys.exit()
#-------------------------END DATA FILTERS--------------------------------------


#### figure ###############################################################
fig = plt.figure()

x = df2['R']
y = df2['Encounters']
z = df2['DormFreq']

ax1 = fig.add_subplot(1, 1, 1, projection = '3d')

ax1.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.0, antialiased=True)
#ax1.set_xlabel('Resources')
ax1.set_ylabel('Encounters')
ax1.set_zlabel('%Dormancy')
plt.tick_params(
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off

plt.show()
#plt.close()
