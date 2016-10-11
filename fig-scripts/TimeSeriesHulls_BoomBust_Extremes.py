from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import statsmodels.api as sm
from random import randint, sample
import os
import sys
from math import isnan, isinf


def plotlines(fig, dat, datlabel, clr, lwidth, sims):

    for sim in sims:
        d = dat[dat['sim'] == sim]
        x = d['ct'].tolist()
        y = np.log10(d[datlabel]).tolist()

        lowess = sm.nonparametric.lowess(y, x, frac=0.01)
        plt.plot(lowess[:, 0], lowess[:, 1], c = clr, ls = '-', lw=lwidth)

    return fig


def plot_hulls(fig, mdat, ylab, clr, clim = 97.5):

    simlist = list(set(mdat['sim'].tolist()))
    print len(simlist)
    sims = sample(simlist, 100)
    X = []
    Y = []

    for i, sim in enumerate(sims):
        d = mdat[mdat['sim'] == sim]
        x = d['ct'].tolist()
        y = np.log10(d[ylab]).tolist()

        for ii, val in enumerate(y):
            if isnan(val) == False and isinf(val) == False:
                X.append(x[ii])
                Y.append(val)

    xran = np.arange(1, 2000, 1).tolist()
    binned = np.digitize(X, xran).tolist()
    bins = [list([]) for _ in xrange(2000)]

    for i, val in enumerate(binned):
        bins[val-1].append(Y[i])

    pct5 = []
    pct95 = []
    xran = []
    for i, _bin in enumerate(bins):
        if len(_bin) > 0:
            pct5.append(np.percentile(_bin, 100 - clim))
            pct95.append(np.percentile(_bin, clim))
            xran.append(i+1)

    plt.fill_between(xran, pct5, pct95, facecolor= clr, alpha=0.5, lw=0.2)
    return fig


mydir = os.path.expanduser('~/GitHub/Micro-Encounter')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")
df = pd.read_csv(mydir + '/results/simulated_data/SimData-TSA-Extremes.csv')

#-------------------------DATA TRANSFORMATIONS -----------------------

df['%Dormant'] = 100 * df['DormantN']/df['N']
df['R'] = df['R'] + 1
df['Encounters'] = df['Encounters'] + 1

#------------------------- DATA FILTERS -------------------

dat = pd.DataFrame(df)

dat1 = dat[dat['MaxRPF'] == 0.001]
dat1 = dat1[dat1['MaxMainFactor'] == 1000]

dat2 = dat[dat['MaxRPF'] == 0.1]
dat2 = dat2[dat2['MaxMainFactor'] == 10]

#-------------------------------------------------------------------------------

#### plot figure ###############################################################
fig = plt.figure()
fs = 8 # fontsize
xlimit = 2000
lwidth = 1
nlim = 1
xlab = 'generations'

#### PLOT 1 #################################################################
fig.add_subplot(2, 2, 1)

plt.plot([0, 1], [-10, -10], c = '0.6', ls = '-', lw=10, label='Dormant')
plt.plot([0, 1], [-10, -10], c = 'm', ls = '-', lw=10, label='Resources')
plt.plot([0, 1], [-10, -10], c = 'c', ls = '-', lw=10, label='$N$')

fig = plot_hulls(fig, dat1, 'DormantN', '0.3')
fig = plot_hulls(fig, dat1, 'R', 'm')
fig = plot_hulls(fig, dat1, 'N', 'c')

plt.xlabel(xlab, fontsize=fs+5)
plt.xlim(0, xlimit)
plt.ylim(0.0, 3.1)
plt.text(150, 2.7, 'Strong dormancy response')
plt.tick_params(axis='both', which='major', labelsize=fs)
plt.legend(bbox_to_anchor=(-0.04, 1.05, 2.48, .2), loc=10, ncol=3, mode="expand",prop={'size':fs+5})


#### PLOT 2 #################################################################
fig.add_subplot(2, 2, 2)

fig = plot_hulls(fig, dat2, 'DormantN', '0.3')
fig = plot_hulls(fig, dat2, 'R', 'm')
fig = plot_hulls(fig, dat2, 'N', 'c')

plt.xlabel(xlab, fontsize=fs+5)
plt.xlim(0, xlimit)
plt.ylim(0.0, 3.1)
plt.text(150, 2.7, 'Weak dormancy response')
plt.tick_params(axis='both', which='major', labelsize=fs)


#### PLOT 3 #################################################################
fig.add_subplot(2, 2, 3)

plt.plot([0, 1], [-10, -10], c = '0.3', ls = '-', lw=2, label='Dormant')
plt.plot([0, 1], [-10, -10], c = 'm', ls = '-', lw=2, label='Resources')
plt.plot([0, 1], [-10, -10], c = 'c', ls = '-', lw=2, label='$N$')

xlab = 'time'
#ylab = '% Dormant'

simlist = list(set(dat1['sim'].tolist()))
print len(simlist)
sims = sample(simlist, nlim)

fig = plotlines(fig, dat1, 'DormantN', '0.3', lwidth, sims)
#fig = plotlines(fig, dat1, 'Encounters', 'm', lwidth, sims)
fig = plotlines(fig, dat1, 'R', 'm', lwidth, sims)
fig = plotlines(fig, dat1, 'N', 'c', lwidth, sims)

#plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.xlim(0, xlimit)
plt.ylim(0.0, 3.1)
plt.text(150, 2.7, 'Strong dormancy response')
plt.tick_params(axis='both', which='major', labelsize=fs)
plt.legend(bbox_to_anchor=(-0.04, 1.05, 2.48, .2), loc=10, ncol=3, mode="expand",prop={'size':fs+5})


#### PLOT 4 #################################################################
fig.add_subplot(2, 2, 4)

xlab = 'time'
#ylab = '$log$'+r'$_{10}$'+'($N$)'

simlist = list(set(dat2['sim'].tolist()))
print len(simlist)
sims = sample(simlist, nlim)

fig = plotlines(fig, dat2, 'DormantN', '0.3', lwidth, sims)
#fig = plotlines(fig, dat2, 'Encounters', 'm', lwidth, sims)
fig = plotlines(fig, dat2, 'R', 'm', lwidth, sims)
fig = plotlines(fig, dat2, 'N', 'c', lwidth, sims)

#plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.xlim(0, xlimit)
plt.ylim(0.0, 3.1)
plt.text(150, 2.7, 'Weak dormancy response')
plt.tick_params(axis='both', which='major', labelsize=fs)


#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.6)
plt.savefig(mydir + '/results/figures/TimeSeriesHulls-BoomBust-Extremes.png', dpi=600, bbox_inches = "tight")
plt.close()
