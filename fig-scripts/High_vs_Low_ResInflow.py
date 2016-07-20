from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys

import scipy as sc
from scipy import stats

mydir = os.path.expanduser('~/GitHub/Micro-Encounter')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")
dat = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')
dat = dat.convert_objects(convert_numeric=True).dropna()


def r_p_bootstrap(x, y, dat):

    rs = []
    ps = []

    for i in range(1):

        rows, cols = dat.shape
        dat = dat.sample(n=rows, replace = False)
        x1 = dat[x]
        y1 = dat[y]

        m, b, r, p, sderr = stats.linregress(x1, y1)
        rs.append(r)
        ps.append(p)

    return [np.mean(rs), np.mean(ps)]


#-------------------------DATA TRANSFORMS---------------------------------------

#dat['DormFreq'] = dat['MeanDormFreq']
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

print dat.shape
#sys.exit()

#-------------------------END DATA TRANSFORMS-----------------------------------

#-------------------------DATA FILTERS------------------------------------------

dat = dat[dat['ResourceComplexityLevel'] != 3]
#dat = dat[dat['TrophicComplexityLevel'] == 3]
dat = dat[dat['SpatialComplexityLevel'] != 1]
dat_high = dat[dat['ResInflow'] > 22] 
dat_low = dat[dat['ResInflow'] < 10] 
#dat = dat[dat['MaxMetMaint'] > 0.002]

print 'size of dat_high:', dat_high.shape
print 'size of dat_low:', dat_low.shape

#-------------------------END DATA FILTERS--------------------------------------




#### figure ###############################################################

fs = 8 # fontsize
fig = plt.figure()
gd = 20
fs = 8
mct = 1

#### PLOT 1 #################################################################
ax = fig.add_subplot(3, 3, 1)

xlab = 'Encounters, '+'$log$'+r'$_{10}$'
ylab = 'Total abundance, '+'$log$'+r'$_{10}$'

x = 'Encounters'
y = 'TotalAbundance'
r, p = r_p_bootstrap(x, y, dat_high)
x = dat_high[x]
y = dat_high[y]

plt.hexbin(x, y, mincnt=mct, gridsize = gd, bins='log', cmap=plt.cm.jet)
plt.ylabel(ylab, fontsize=fs+1)
plt.xlabel(xlab, fontsize=fs+1)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
#plt.ylim(0.5, 3.5)
plt.text(0.8, 3.4, 'r = '+str(round(r,2)), fontsize=fs+1, color='k')

plt.text(-0.4, 3.25, 'High Inflow', rotation=90, fontsize=fs+6, color='k')
plt.text(-0.4, 1.0, 'Low Inflow', rotation=90, fontsize=fs+6, color='k')

#### PLOT 2 #################################################################
ax = fig.add_subplot(3, 3, 2)

xlab = 'Encounters, '+'$log$'+r'$_{10}$'
ylab = 'Productivity, '+'$log$'+r'$_{10}$'

x = 'Encounters'
y = 'Production'
r, p = r_p_bootstrap(x, y, dat_high)
x = dat_high[x]
y = dat_high[y]

plt.hexbin(x, y, mincnt=mct, gridsize = gd, bins='log', cmap=plt.cm.jet)
plt.ylabel(ylab, fontsize=fs+1)
plt.xlabel(xlab, fontsize=fs+1)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
#plt.ylim(-1.0, 2.0)
plt.text(1.0, 1.3, 'r = '+str(round(r,2)), fontsize=fs+1, color='k')

#### PLOT 3 #################################################################
ax = fig.add_subplot(3, 3, 3)

xlab = 'Encounters, '+'$log$'+r'$_{10}$'
ylab = '% Dormancy, '+'$log$'+r'$_{10}$'

x = 'Encounters'
y = 'DormFreq'
r, p = r_p_bootstrap(x, y, dat_high)
x = dat_high[x]
y = dat_high[y]

plt.hexbin(x, y, mincnt=mct, gridsize = gd, bins='log', cmap=plt.cm.jet)
plt.ylabel(ylab, fontsize=fs+1)
plt.xlabel(xlab, fontsize=fs+1)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.xlim(0.5, 3.0)
#plt.ylim(0.1, 1.2)
plt.text(0.7, -0.72, 'r = '+str(round(r,2)), fontsize=fs+1, color='k')

#### PLOT 4 #################################################################
ax = fig.add_subplot(3, 3, 4)

xlab = 'Encounters, '+'$log$'+r'$_{10}$'
ylab = 'Total abundance, '+'$log$'+r'$_{10}$'

x = 'Encounters'
y = 'TotalAbundance'
r, p = r_p_bootstrap(x, y, dat_low)
x = dat_low[x]
y = dat_low[y]

plt.hexbin(x, y, mincnt=mct, gridsize = gd, bins='log', cmap=plt.cm.jet)
plt.ylabel(ylab, fontsize=fs+1)
plt.xlabel(xlab, fontsize=fs+1)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
#plt.ylim(0.5, 3.5)
plt.text(-0.55, 2.5, 'r = '+str(round(r,2)), fontsize=fs+1, color='k')

#### PLOT 5 #################################################################
ax = fig.add_subplot(3, 3, 5)

xlab = 'Encounters, '+'$log$'+r'$_{10}$'
ylab = 'Productivity, '+'$log$'+r'$_{10}$'

x = 'Encounters'
y = 'Production'
r, p = r_p_bootstrap(x, y, dat_low)
x = dat_low[x]
y = dat_low[y]

plt.hexbin(x, y, mincnt=mct, gridsize = gd, bins='log', cmap=plt.cm.jet)
plt.ylabel(ylab, fontsize=fs+1)
plt.xlabel(xlab, fontsize=fs+1)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
#plt.ylim(-1.0, 2.0)
plt.text(-0.5, 0.25, 'r = '+str(round(r,2)), fontsize=fs+1, color='k')

#### PLOT 6 #################################################################
ax = fig.add_subplot(3, 3, 6)

xlab = 'Encounters, '+'$log$'+r'$_{10}$'
ylab = '% Dormancy, '+'$log$'+r'$_{10}$'

x = 'Encounters'
y = 'DormFreq'
r, p = r_p_bootstrap(x, y, dat_low)
x = dat_low[x]
y = dat_low[y]

plt.hexbin(x, y, mincnt=mct, gridsize = gd, bins='log', cmap=plt.cm.jet)
plt.ylabel(ylab, fontsize=fs+1)
plt.xlabel(xlab, fontsize=fs+1)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.ylim(-1.2, 0.0)
plt.xlim(-1.1, 2.1)
plt.text(-0.95, -1.0, 'r = '+str(round(r,2)), fontsize=fs+1, color='k')


#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.6, hspace=0.5)
plt.savefig(mydir + '/results/figures/High_vs_Low_ResInflow.png', dpi=600, bbox_inches = "tight")
#plt.show()
#plt.close()
