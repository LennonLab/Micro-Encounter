from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from random import randint
import os
import sys
import statsmodels.api as sm
from math import isnan, isinf


def randcolor():
    c1 = randint(0,255)
    c2 = randint(0,255)
    c3 = randint(0,255)

    clr = '#%02x%02x%02x' % (c1, c2, c3)
    return clr

def plot_dat(fig, mdat, ylab, clr, clim = 65):

    simlist = list(set(mdat['sim'].tolist()))
    print len(simlist)
    X = []
    Y = []

    for i, sim in enumerate(simlist):
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
df = pd.read_csv(mydir + '/results/simulated_data/SimData-TSA.csv')

#-------------------------DATA TRANSFORMATIONS -----------------------

df['%Dormant'] = 100*df['DormantN']/df['N']

#------------------------- DATA FILTERS -------------------

dat = pd.DataFrame(df)
dat = dat[~dat['SpatialComplexity'].str.contains('-wellmixed-')]
#dat = dat[dat['ct'] > 500]
#dat = dat[dat['Immigration'] < 0.05]
#dat = dat[dat['MaxMaint'] >= 5]
#dat = dat[dat['MaxMainFactor'] < 30]
#dat = dat[dat['width'] <= 250]

dat = dat[~dat['ResourceComplexity'].str.contains('-simple-')]
dat = dat[~dat['ResourceComplexity'].str.contains('-monoculture-')]
dat = dat[~dat['TrophicComplexity'].str.contains('-scavenging-')]
dat = dat[~dat['TrophicComplexity'].str.contains('-crossfeeding-')]

dat1 = dat[dat['SpatialComplexity'].str.contains('-chemotaxis-')]
dat2 = dat[dat['SpatialComplexity'].str.contains('-randwalk-')]
dat3 = dat[dat['SpatialComplexity'].str.contains('-none-')]
dat4 = dat[dat['SpatialComplexity'].str.contains('-none-')]

#-------------------------------------------------------------------------------

#### plot figure ###############################################################
fig = plt.figure()
fs = 8 # fontsize
xlimit = 2000
xlab = 'Generations'

#### PLOT 1 #################################################################
fig.add_subplot(2, 2, 1)

plt.plot([0, 1], [-10, -10], c = '0.3', ls = '-', alpha = 0.5, lw=10, label='Seed Bank')
#plt.plot([0, 1], [-10, -10], c = 'm', ls = '-', alpha = 0.5, lw=10, label='Resources')
plt.plot([0, 1], [-10, -10], c = 'c', ls = '-', alpha = 0.5, lw=10, label='Encounters')

ylab = 'Number'
fig = plot_dat(fig, dat1, '%Dormant', '0.3')
fig = plot_dat(fig, dat1, 'Encounters', 'c')

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.ylim(0.2, 2)
plt.tick_params(axis='both', which='major', labelsize=fs)
plt.legend(bbox_to_anchor=(-0.04, 1.05, 2.48, .2), loc=10, ncol=3, mode="expand",prop={'size':fs+5})


#### PLOT 2 #################################################################
fig.add_subplot(2, 2, 2)

datylab = 'Number'
fig = plot_dat(fig, dat2, '%Dormant', '0.3')
fig = plot_dat(fig, dat2, 'Encounters', 'c')

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
#plt.ylim(0.0, 3.1)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### PLOT 3 #################################################################
fig.add_subplot(2, 2, 3)


datylab = 'Number'
fig = plot_dat(fig, dat3, '%Dormant', '0.3')
fig = plot_dat(fig, dat3, 'Encounters', 'c')

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
#plt.ylim(-0.1, 3.0)
plt.tick_params(axis='both', which='major', labelsize=fs)


#### PLOT 4 #################################################################
fig.add_subplot(2, 2, 4)


datylab = 'Number'
fig = plot_dat(fig, dat3, '%Dormant', '0.3')
fig = plot_dat(fig, dat3, 'R', 'c')

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
#plt.ylim(-0.1, 2.1)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/TimeSeries-BoomBust.png', dpi=600, bbox_inches = "tight")
plt.close()
