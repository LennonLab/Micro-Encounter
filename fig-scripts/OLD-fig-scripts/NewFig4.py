from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys

import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import summary_table

mydir = os.path.expanduser('~/GitHub/Micro-Encounter')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")
dat = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')
dat = dat.convert_objects(convert_numeric=True).dropna()

#print dat.shape
#sys.exit()

#-------------------------DATA FILTERS------------------------------------------

dat = dat[dat['SpatialComplexityLevel'] == 1]
dat = dat[dat['ResourceComplexityLevel'] <= 1]
#dat = dat[dat['TrophicComplexityLevel'] <= 3]

dat['DormFreq'] = np.log10(dat['MeanDormFreq'])
dat = dat[np.isfinite(dat['DormFreq'])]

dat['Encounters'] = np.log10(dat['MeanEncounter'])
dat = dat[np.isfinite(dat['Encounters'])]

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
dat1 = dat[dat['TrophicComplexityLevel'] == 1]
dat2 = dat[dat['TrophicComplexityLevel'] == 2]
dat3 = dat[dat['TrophicComplexityLevel'] == 3]

label1 = 'No trophic complexity'
label2 = 'Crossfeeding'
label3 = 'Recycling'

#### PLOT 1 #################################################################
fig.add_subplot(2, 2, 1)
xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
ylab = '% Dormancy, '+'$log$'+r'$_{10}$'
width = 1

f = smf.ols('DormFreq ~ Encounters', dat1).fit()
st, data, ss2 = summary_table(f, alpha=0.05)
mean_ci_low, mean_ci_upp = data[:,4:6].T
x = dat1['Encounters'].tolist()
y = dat1['DormFreq'].tolist()
ylow = mean_ci_low.tolist()
yupp = mean_ci_upp.tolist()
x, y, yupp, ylow = zip(*sorted(zip(x, y, yupp, ylow)))
plt.scatter(x, y, color = 'm', alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = 'm', alpha = 0.9, lw=width, label=label1)
plt.plot(x, yupp, color = 'm', alpha = 0.9, lw=width)


f = smf.ols('DormFreq ~ Encounters', dat2).fit()
st, data, ss2 = summary_table(f, alpha=0.05)
mean_ci_low, mean_ci_upp = data[:,4:6].T
x = dat2['Encounters'].tolist()
y = dat2['DormFreq'].tolist()
ylow = mean_ci_low.tolist()
yupp = mean_ci_upp.tolist()
x, y, yupp, ylow = zip(*sorted(zip(x, y, yupp, ylow)))
plt.scatter(x, y, color = 'steelblue', alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = 'steelblue', alpha = 0.9, lw=width, label=label2)
plt.plot(x, yupp, color = 'steelblue', alpha = 0.9, lw=width)


f = smf.ols('DormFreq ~ Encounters', dat3).fit()
st, data, ss2 = summary_table(f, alpha=0.05)
mean_ci_low, mean_ci_upp = data[:,4:6].T
x = dat3['Encounters'].tolist()
y = dat3['DormFreq'].tolist()
ylow = mean_ci_low.tolist()
yupp = mean_ci_upp.tolist()
x, y, yupp, ylow = zip(*sorted(zip(x, y, yupp, ylow)))
plt.scatter(x, y, color = 'goldenrod', alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = 'goldenrod', alpha = 0.9, lw=width, label=label3)
plt.plot(x, yupp, color = 'goldenrod', alpha = 0.9, lw=width)


plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.ylim(-1.0, 0.1)

plt.tick_params(axis='both', which='major', labelsize=fs)
plt.legend(bbox_to_anchor=(-0.04, 1.05, 2.48, .2), loc=10, ncol=3, mode="expand",prop={'size':fs})

#### PLOT 2 ################################
fig.add_subplot(2, 2, 2)

xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
ylab = 'Individual production, '+'$log$'+r'$_{10}$'

f = smf.ols('Production ~ Encounters', dat1).fit()
st, data, ss2 = summary_table(f, alpha=0.05)
mean_ci_low, mean_ci_upp = data[:,4:6].T
x = dat1['Encounters'].tolist()
y = dat1['Production'].tolist()
ylow = mean_ci_low.tolist()
yupp = mean_ci_upp.tolist()
x, y, yupp, ylow = zip(*sorted(zip(x, y, yupp, ylow)))
plt.scatter(x, y, color = 'm', alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = 'm', alpha = 0.9, lw=width, label=label1)
plt.plot(x, yupp, color = 'm', alpha = 0.9, lw=width)


f = smf.ols('Production ~ Encounters', dat2).fit()
st, data, ss2 = summary_table(f, alpha=0.05)
mean_ci_low, mean_ci_upp = data[:,4:6].T
x = dat2['Encounters'].tolist()
y = dat2['Production'].tolist()
ylow = mean_ci_low.tolist()
yupp = mean_ci_upp.tolist()
x, y, yupp, ylow = zip(*sorted(zip(x, y, yupp, ylow)))
plt.scatter(x, y, color = 'steelblue', alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = 'steelblue', alpha = 0.9, lw=width, label=label2)
plt.plot(x, yupp, color = 'steelblue', alpha = 0.9, lw=width)


f = smf.ols('Production ~ Encounters', dat3).fit()
st, data, ss2 = summary_table(f, alpha=0.05)
mean_ci_low, mean_ci_upp = data[:,4:6].T
x = dat3['Encounters'].tolist()
y = dat3['Production'].tolist()
ylow = mean_ci_low.tolist()
yupp = mean_ci_upp.tolist()
x, y, yupp, ylow = zip(*sorted(zip(x, y, yupp, ylow)))
plt.scatter(x, y, color = 'goldenrod', alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = 'goldenrod', alpha = 0.9, lw=width, label=label3)
plt.plot(x, yupp, color = 'goldenrod', alpha = 0.9, lw=width)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.ylim(-2.0, 2.0)
#plt.xlim(0.1, 300)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### PLOT 3 #################################################################
fig.add_subplot(2, 2, 3)

xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
ylab = 'Total abundance, '+'$log$'+r'$_{10}$'

f = smf.ols('TotalAbundance ~ Encounters', dat1).fit()
st, data, ss2 = summary_table(f, alpha=0.05)
mean_ci_low, mean_ci_upp = data[:,4:6].T
x = dat1['Encounters'].tolist()
y = dat1['TotalAbundance'].tolist()
ylow = mean_ci_low.tolist()
yupp = mean_ci_upp.tolist()
x, y, yupp, ylow = zip(*sorted(zip(x, y, yupp, ylow)))
plt.scatter(x, y, color = 'm', alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = 'm', alpha = 0.9, lw=width, label=label1)
plt.plot(x, yupp, color = 'm', alpha = 0.9, lw=width)


f = smf.ols('TotalAbundance ~ Encounters', dat2).fit()
st, data, ss2 = summary_table(f, alpha=0.05)
mean_ci_low, mean_ci_upp = data[:,4:6].T
x = dat2['Encounters'].tolist()
y = dat2['TotalAbundance'].tolist()
ylow = mean_ci_low.tolist()
yupp = mean_ci_upp.tolist()
x, y, yupp, ylow = zip(*sorted(zip(x, y, yupp, ylow)))
plt.scatter(x, y, color = 'steelblue', alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = 'steelblue', alpha = 0.9, lw=width, label=label2)
plt.plot(x, yupp, color = 'steelblue', alpha = 0.9, lw=width)


f = smf.ols('TotalAbundance ~ Encounters', dat3).fit()
st, data, ss2 = summary_table(f, alpha=0.05)
mean_ci_low, mean_ci_upp = data[:,4:6].T
x = dat3['Encounters'].tolist()
y = dat3['TotalAbundance'].tolist()
ylow = mean_ci_low.tolist()
yupp = mean_ci_upp.tolist()
x, y, yupp, ylow = zip(*sorted(zip(x, y, yupp, ylow)))
plt.scatter(x, y, color = 'goldenrod', alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = 'goldenrod', alpha = 0.9, lw=width, label=label3)
plt.plot(x, yupp, color = 'goldenrod', alpha = 0.9, lw=width)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
#plt.xlim(0.15, 300)
plt.ylim(0.5, 3.1)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### PLOT 4 #################################################################
fig.add_subplot(2, 2, 4)

xlab = 'Average encounters, '+'$log$'+r'$_{10}$'
ylab = 'Active abundance, '+'$log$'+r'$_{10}$'


f = smf.ols('ActiveAbundance ~ Encounters', dat1).fit()
st, data, ss2 = summary_table(f, alpha=0.05)
mean_ci_low, mean_ci_upp = data[:,4:6].T
x = dat1['Encounters'].tolist()
y = dat1['ActiveAbundance'].tolist()
ylow = mean_ci_low.tolist()
yupp = mean_ci_upp.tolist()
x, y, yupp, ylow = zip(*sorted(zip(x, y, yupp, ylow)))
plt.scatter(x, y, color = 'm', alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = 'm', alpha = 0.9, lw=width, label=label1)
plt.plot(x, yupp, color = 'm', alpha = 0.9, lw=width)


f = smf.ols('ActiveAbundance ~ Encounters', dat2).fit()
st, data, ss2 = summary_table(f, alpha=0.05)
mean_ci_low, mean_ci_upp = data[:,4:6].T
x = dat2['Encounters'].tolist()
y = dat2['ActiveAbundance'].tolist()
ylow = mean_ci_low.tolist()
yupp = mean_ci_upp.tolist()
x, y, yupp, ylow = zip(*sorted(zip(x, y, yupp, ylow)))
plt.scatter(x, y, color = 'steelblue', alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = 'steelblue', alpha = 0.9, lw=width, label=label2)
plt.plot(x, yupp, color = 'steelblue', alpha = 0.9, lw=width)


f = smf.ols('ActiveAbundance ~ Encounters', dat3).fit()
st, data, ss2 = summary_table(f, alpha=0.05)
mean_ci_low, mean_ci_upp = data[:,4:6].T
x = dat3['Encounters'].tolist()
y = dat3['ActiveAbundance'].tolist()
ylow = mean_ci_low.tolist()
yupp = mean_ci_upp.tolist()
x, y, yupp, ylow = zip(*sorted(zip(x, y, yupp, ylow)))
plt.scatter(x, y, color = 'goldenrod', alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = 'goldenrod', alpha = 0.9, lw=width, label=label3)
plt.plot(x, yupp, color = 'goldenrod', alpha = 0.9, lw=width)


plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
#plt.xlim(0.15, 1000)
plt.ylim(-0.5, 3.1)
plt.tick_params(axis='both', which='major', labelsize=fs)


plt.savefig(mydir + '/results/figures/Fig4-Trophic_TC1&2-RC3-SC1-.png', dpi=600, bbox_inches = "tight")
#plt.show()
plt.close()
