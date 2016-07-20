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

    for i in range(100):

        dat = dat.sample(n=100, replace = True)
        x1 = dat[x]
        y1 = dat[y]

        m, b, r, p, sderr = stats.linregress(x1, y1)
        rs.append(r)
        ps.append(p)

    return [np.mean(rs), np.mean(ps)]


#-------------------------DATA FILTERS------------------------------------------

#dat = dat[dat['ResourceComplexityLevel'] != 3]
#dat = dat[dat['TrophicComplexityLevel'] == 3]
#dat = dat[dat['SpatialComplexityLevel'] != 1]
#dat = dat[dat['MeanTotalAbundance'] >= 20]
#dat = dat[dat['MaxMetMaint'] > 0.002]

dat = dat[dat['height'] <= 6]
dat = dat[dat['width'] <= 6]


#-------------------------END DATA FILTERS--------------------------------------


#-------------------------DATA TRANSFORMS---------------------------------------

dat['DormFreq'] = np.log10(dat['MeanDormFreq'])
dat = dat[np.isfinite(dat['DormFreq'])]

dat['Encounters'] = np.log10(dat['MeanEncounter'])
dat = dat[np.isfinite(dat['Encounters'])]

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


#### figure ###############################################################

fs = 8 # fontsize
fig = plt.figure()
gd = 20
fs = 8
mct = 1

#### PLOT 1 #################################################################
ax = fig.add_subplot(3, 3, 1)

xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
ylab = 'Total abundance, '+'$log$'+r'$_{10}$'

x = 'Encounters'
y = 'TotalAbundance'
r, p = r_p_bootstrap(x, y, dat)
x = dat[x]
y = dat[y]

plt.hexbin(x, y, mincnt=mct, gridsize = gd, bins='log', cmap=plt.cm.jet)
plt.ylabel(ylab, fontsize=fs+1)
plt.xlabel(xlab, fontsize=fs+1)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
#plt.ylim(0.5, 2.5)
plt.text(-1.0, 2.1, 'r = '+str(round(r,2)), fontsize=fs+1, color='k')

#### PLOT 2 #################################################################
ax = fig.add_subplot(3, 3, 2)

xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
ylab = 'Productivity, '+'$log$'+r'$_{10}$'

x = 'Encounters'
y = 'Production'
r, p = r_p_bootstrap(x, y, dat)
x = dat[x]
y = dat[y]

plt.hexbin(x, y, mincnt=mct, gridsize = gd, bins='log', cmap=plt.cm.jet)
plt.ylabel(ylab, fontsize=fs+1)
plt.xlabel(xlab, fontsize=fs+1)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
#plt.ylim(-1.0, 1.5)
plt.text(-1.0, 1.3, 'r = '+str(round(r,2)), fontsize=fs+1, color='k')

#### PLOT 3 #################################################################
ax = fig.add_subplot(3, 3, 3)

xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
ylab = '% Dormancy, '+'$log$'+r'$_{10}$'

x = 'Encounters'
y = 'DormFreq'
r, p = r_p_bootstrap(x, y, dat)
x = dat[x]
y = dat[y]

plt.hexbin(x, y, mincnt=mct, gridsize = gd, bins='log', cmap=plt.cm.jet)
plt.ylabel(ylab, fontsize=fs+1)
plt.xlabel(xlab, fontsize=fs+1)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
#plt.ylim(0.1,1.2)
plt.text(-1.0, 1.05, 'r = '+str(round(r,2)), fontsize=fs+1, color='k')

#### PLOT 4 #################################################################
ax = fig.add_subplot(3, 3, 4)

xlab = 'Resource inflow, '+'$log$'+r'$_{10}$'
ylab = 'Total abundance, '+'$log$'+r'$_{10}$'

x = 'Res_Inflow'
y = 'TotalAbundance'
r, p = r_p_bootstrap(x, y, dat)
x = dat[x]
y = dat[y]

plt.hexbin(x, y, mincnt=mct, gridsize = gd, bins='log', cmap=plt.cm.jet)
plt.ylabel(ylab, fontsize=fs+1)
plt.xlabel(xlab, fontsize=fs+1)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
#plt.ylim(0.5,3.5)
plt.text(0.5, 3.1, 'r = '+str(round(r,2)), fontsize=fs+1, color='k')

#### PLOT 5 #################################################################
ax = fig.add_subplot(3, 3, 5)

xlab = 'Resource inflow, '+'$log$'+r'$_{10}$'
ylab = 'Productivity, '+'$log$'+r'$_{10}$'

x = 'Res_Inflow'
y = 'Production'
r, p = r_p_bootstrap(x, y, dat)
x = dat[x]
y = dat[y]

plt.hexbin(x, y, mincnt=mct, gridsize = gd, bins='log', cmap=plt.cm.jet)
plt.ylabel(ylab, fontsize=fs+1)
plt.xlabel(xlab, fontsize=fs+1)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
#plt.ylim(-1, 2)
plt.text(0.1, 1.6, 'r = '+str(round(r,2)), fontsize=fs+1, color='k')

#### PLOT 6 #################################################################
ax = fig.add_subplot(3, 3, 6)

xlab = 'Resource inflow, '+'$log$'+r'$_{10}$'
ylab = '% Dormancy, '+'$log$'+r'$_{10}$'

x = 'Res_Inflow'
y = 'DormFreq'
r, p = r_p_bootstrap(x, y, dat)
x = dat[x]
y = dat[y]

plt.hexbin(x, y, mincnt=mct, gridsize = gd, bins='log', cmap=plt.cm.jet)
plt.ylabel(ylab, fontsize=fs+1)
plt.xlabel(xlab, fontsize=fs+1)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
#plt.ylim(0.1,1.2)
plt.text(0.1, 1.05, 'r = '+str(round(r,2)), fontsize=fs+1, color='k')

#### PLOT 7 #################################################################
ax = fig.add_subplot(3, 3, 7)

xlab = 'Total resources , '+'$log$'+r'$_{10}$'
ylab = 'Total abundance, '+'$log$'+r'$_{10}$'

x = 'AvgResourceParticles'
y = 'TotalAbundance'
r, p = r_p_bootstrap(x, y, dat)
x = dat[x]
y = dat[y]

plt.hexbin(x, y, mincnt=mct, gridsize = gd, bins='log', cmap=plt.cm.jet)
plt.ylabel(ylab, fontsize=fs+1)
plt.xlabel(xlab, fontsize=fs+1)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
#plt.ylim(0.5,3.5)
plt.text(0.5, 3.1, 'r = '+str(round(r,2)), fontsize=fs+1, color='k')

#### PLOT 8 #################################################################
ax = fig.add_subplot(3, 3, 8)

xlab = 'Total resources , '+'$log$'+r'$_{10}$'
ylab = 'Productivity, '+'$log$'+r'$_{10}$'

x = 'AvgResourceParticles'
y = 'Production'
r, p = r_p_bootstrap(x, y, dat)
x = dat[x]
y = dat[y]

plt.hexbin(x, y, mincnt=mct, gridsize = gd, bins='log', cmap=plt.cm.jet)
plt.ylabel(ylab, fontsize=fs+1)
plt.xlabel(xlab, fontsize=fs+1)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
#plt.ylim(-1, 2)
plt.text(0.5, 1.6, 'r = '+str(round(r,2)), fontsize=fs+1, color='k')

#### PLOT 9 #################################################################
ax = fig.add_subplot(3, 3, 9)

xlab = 'Total resources , '+'$log$'+r'$_{10}$'
ylab = '% Dormancy, '+'$log$'+r'$_{10}$'

x = 'AvgResourceParticles'
y = 'DormFreq'
r, p = r_p_bootstrap(x, y, dat)
x = dat[x]
y = dat[y]

plt.hexbin(x, y, mincnt=mct, gridsize = gd, bins='log', cmap=plt.cm.jet)
plt.ylabel(ylab, fontsize=fs+1)
plt.xlabel(xlab, fontsize=fs+1)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
#plt.ylim(0.1,1.2)
plt.text(0.5, 1.05, 'r = '+str(round(r,2)), fontsize=fs+1, color='k')

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.6, hspace=0.5)
#plt.savefig(mydir + '/results/figures/Res_3by3.png', dpi=600, bbox_inches = "tight")
plt.show()
#plt.close()
