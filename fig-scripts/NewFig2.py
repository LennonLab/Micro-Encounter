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
dat = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')
dat = dat.convert_objects(convert_numeric=True).dropna()

color1 = 'm'
color2 = 'steelblue'
color3 = 'goldenrod'

#-------------------------DATA FILTERS------------------------------------------

dat = dat[dat['ResourceComplexityLevel'] == 2]
dat = dat[dat['TrophicComplexityLevel'] == 1]

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

dat1 = dat[dat['SpatialComplexityLevel'] == 1]
dat2 = dat[dat['SpatialComplexityLevel'] == 2]
dat3 = dat[dat['SpatialComplexityLevel'] == 3]

label1 = 'White noise'
label2 = 'Aggregated w/ Random walks'
label3 = 'Aggregated w/ chemotaxis'

Xs = []
Ys = []
colors = []

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
Xs.extend(x)
Ys.extend(y)
colors.extend([color1]*len(x))
plt.plot(x, ylow, color = color1, alpha = 0.9, lw=width, label=label1)
plt.plot(x, yupp, color = color1, alpha = 0.9, lw=width)


f = smf.ols('DormFreq ~ Encounters', dat2).fit()
st, data, ss2 = summary_table(f, alpha=0.05)
mean_ci_low, mean_ci_upp = data[:,4:6].T
x = dat2['Encounters'].tolist()
y = dat2['DormFreq'].tolist()
ylow = mean_ci_low.tolist()
yupp = mean_ci_upp.tolist()
x, y, yupp, ylow = zip(*sorted(zip(x, y, yupp, ylow)))
Xs.extend(x)
Ys.extend(y)
colors.extend([color2]*len(x))
plt.plot(x, ylow, color = color2, alpha = 0.9, lw=width, label=label2)
plt.plot(x, yupp, color = color2, alpha = 0.9, lw=width)


f = smf.ols('DormFreq ~ Encounters', dat3).fit()
st, data, ss2 = summary_table(f, alpha=0.05)
mean_ci_low, mean_ci_upp = data[:,4:6].T
x = dat3['Encounters'].tolist()
y = dat3['DormFreq'].tolist()
ylow = mean_ci_low.tolist()
yupp = mean_ci_upp.tolist()
x, y, yupp, ylow = zip(*sorted(zip(x, y, yupp, ylow)))
Xs.extend(x)
Ys.extend(y)
colors.extend([color3]*len(x))
plt.plot(x, ylow, color = color3, alpha = 0.9, lw=width, label=label3)
plt.plot(x, yupp, color = color3, alpha = 0.9, lw=width)

indices = range(0, len(Xs))
shuffle(indices)

for i in indices:
    plt.scatter(Xs[i], Ys[i], color = colors[i], alpha = 0.6 , s = 5, linewidths = 0.0)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.ylim(-1.0, 0.1)
plt.tick_params(axis='both', which='major', labelsize=fs)
plt.legend(bbox_to_anchor=(-0.04, 1.05, 2.48, .2), loc=10, ncol=3, mode="expand",prop={'size':fs})

#### PLOT 2 ################################
fig.add_subplot(2, 2, 2)

Xs = []
Ys = []
colors = []

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
Xs.extend(x)
Ys.extend(y)
colors.extend([color1]*len(x))
plt.scatter(x, y, color = color1, alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = color1, alpha = 0.9, lw=width, label=label1)
plt.plot(x, yupp, color = color1, alpha = 0.9, lw=width)


f = smf.ols('Production ~ Encounters', dat2).fit()
st, data, ss2 = summary_table(f, alpha=0.05)
mean_ci_low, mean_ci_upp = data[:,4:6].T
x = dat2['Encounters'].tolist()
y = dat2['Production'].tolist()
ylow = mean_ci_low.tolist()
yupp = mean_ci_upp.tolist()
x, y, yupp, ylow = zip(*sorted(zip(x, y, yupp, ylow)))
Xs.extend(x)
Ys.extend(y)
colors.extend([color2]*len(x))
plt.scatter(x, y, color = color2, alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = color2, alpha = 0.9, lw=width, label=label2)
plt.plot(x, yupp, color = color2, alpha = 0.9, lw=width)


f = smf.ols('Production ~ Encounters', dat3).fit()
st, data, ss2 = summary_table(f, alpha=0.05)
mean_ci_low, mean_ci_upp = data[:,4:6].T
x = dat3['Encounters'].tolist()
y = dat3['Production'].tolist()
ylow = mean_ci_low.tolist()
yupp = mean_ci_upp.tolist()
x, y, yupp, ylow = zip(*sorted(zip(x, y, yupp, ylow)))
Xs.extend(x)
Ys.extend(y)
colors.extend([color3]*len(x))
plt.scatter(x, y, color = color3, alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = color3, alpha = 0.9, lw=width, label=label3)
plt.plot(x, yupp, color = color3, alpha = 0.9, lw=width)


indices = range(0, len(Xs))
shuffle(indices)

for i in indices:
    plt.scatter(Xs[i], Ys[i], color = colors[i], alpha = 0.6 , s = 5, linewidths = 0.0)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.ylim(-2.0, 2.0)
#plt.xlim(0.1, 300)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### PLOT 3 #################################################################
fig.add_subplot(2, 2, 3)

Xs = []
Ys = []
colors = []

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
Xs.extend(x)
Ys.extend(y)
colors.extend([color1]*len(x))
plt.scatter(x, y, color = color1, alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = color1, alpha = 0.9, lw=width, label=label1)
plt.plot(x, yupp, color = color1, alpha = 0.9, lw=width)


f = smf.ols('TotalAbundance ~ Encounters', dat2).fit()
st, data, ss2 = summary_table(f, alpha=0.05)
mean_ci_low, mean_ci_upp = data[:,4:6].T
x = dat2['Encounters'].tolist()
y = dat2['TotalAbundance'].tolist()
ylow = mean_ci_low.tolist()
yupp = mean_ci_upp.tolist()
x, y, yupp, ylow = zip(*sorted(zip(x, y, yupp, ylow)))
Xs.extend(x)
Ys.extend(y)
colors.extend([color2]*len(x))
plt.scatter(x, y, color = color2, alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = color2, alpha = 0.9, lw=width, label=label2)
plt.plot(x, yupp, color = color2, alpha = 0.9, lw=width)


f = smf.ols('TotalAbundance ~ Encounters', dat3).fit()
st, data, ss2 = summary_table(f, alpha=0.05)
mean_ci_low, mean_ci_upp = data[:,4:6].T
x = dat3['Encounters'].tolist()
y = dat3['TotalAbundance'].tolist()
ylow = mean_ci_low.tolist()
yupp = mean_ci_upp.tolist()
x, y, yupp, ylow = zip(*sorted(zip(x, y, yupp, ylow)))
Xs.extend(x)
Ys.extend(y)
colors.extend([color3]*len(x))
plt.scatter(x, y, color = color3, alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = color3, alpha = 0.9, lw=width, label=label3)
plt.plot(x, yupp, color = color3, alpha = 0.9, lw=width)


indices = range(0, len(Xs))
shuffle(indices)

for i in indices:
    plt.scatter(Xs[i], Ys[i], color = colors[i], alpha = 0.6 , s = 5, linewidths = 0.0)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
#plt.xlim(0.15, 300)
plt.ylim(0.5, 3.1)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### PLOT 4 #################################################################
fig.add_subplot(2, 2, 4)

Xs = []
Ys = []
colors = []

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
Xs.extend(x)
Ys.extend(y)
colors.extend([color1]*len(x))
plt.scatter(x, y, color = color1, alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = color1, alpha = 0.9, lw=width, label=label1)
plt.plot(x, yupp, color = color1, alpha = 0.9, lw=width)


f = smf.ols('ActiveAbundance ~ Encounters', dat2).fit()
st, data, ss2 = summary_table(f, alpha=0.05)
mean_ci_low, mean_ci_upp = data[:,4:6].T
x = dat2['Encounters'].tolist()
y = dat2['ActiveAbundance'].tolist()
ylow = mean_ci_low.tolist()
yupp = mean_ci_upp.tolist()
x, y, yupp, ylow = zip(*sorted(zip(x, y, yupp, ylow)))
Xs.extend(x)
Ys.extend(y)
colors.extend([color2]*len(x))
plt.scatter(x, y, color = color2, alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = color2, alpha = 0.9, lw=width, label=label2)
plt.plot(x, yupp, color = color2, alpha = 0.9, lw=width)


f = smf.ols('ActiveAbundance ~ Encounters', dat3).fit()
st, data, ss2 = summary_table(f, alpha=0.05)
mean_ci_low, mean_ci_upp = data[:,4:6].T
x = dat3['Encounters'].tolist()
y = dat3['ActiveAbundance'].tolist()
ylow = mean_ci_low.tolist()
yupp = mean_ci_upp.tolist()
x, y, yupp, ylow = zip(*sorted(zip(x, y, yupp, ylow)))
Xs.extend(x)
Ys.extend(y)
colors.extend([color3]*len(x))
plt.scatter(x, y, color = color3, alpha = 0.3 , s = 5, linewidths = 0.0)
plt.plot(x, ylow, color = color3, alpha = 0.9, lw=width, label=label3)
plt.plot(x, yupp, color = color3, alpha = 0.9, lw=width)

indices = range(0, len(Xs))
shuffle(indices)

for i in indices:
    plt.scatter(Xs[i], Ys[i], color = colors[i], alpha = 0.6 , s = 5, linewidths = 0.0)


plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
#plt.xlim(0.15, 1000)
plt.ylim(-0.5, 3.1)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/Fig2-Spatial_RC2-TC1.png', dpi=600, bbox_inches = "tight")
#plt.show()
