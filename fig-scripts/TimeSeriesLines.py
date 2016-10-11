from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from random import randint
import os
import sys
import statsmodels.api as sm


def randcolor():
    c1 = randint(0,255)
    c2 = randint(0,255)
    c3 = randint(0,255)

    clr = '#%02x%02x%02x' % (c1, c2, c3)
    return clr


def plotlines(fig, dat, datlabel, clr, lwidth, nlim):

    simlist = list(set(dat['sim'].tolist()))
    print len(simlist)
    for i, sim in enumerate(simlist):
        d = dat[dat['sim'] == sim]
        x = d['ct'].tolist()
        if datlabel == '%Dormant':
            y = d[datlabel].tolist()
        else:
            y = np.log10(d[datlabel]).tolist()

        lowess = sm.nonparametric.lowess(y, x, frac=0.1)
        plt.plot(lowess[:, 0], lowess[:, 1], c = clr, ls = '-', lw=lwidth)
        if i > nlim: break

    return fig

mydir = os.path.expanduser('~/GitHub/Micro-Encounter')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")
df = pd.read_csv(mydir + '/results/simulated_data/SimData-TSA.csv')
#df = df.convert_objects(convert_numeric=True).dropna()

#-------------------------DATA TRANSFORMATIONS -----------------------

df['%Dormant'] = 100 * df['DormantN']/df['N']

#------------------------- DATA FILTERS -------------------

dat = pd.DataFrame(df)
#dat = dat[dat['ct'] > 500]
dat = dat[~dat['SpatialComplexity'].str.contains('-wellmixed-')]
#dat = dat[dat['ResourceComplexity'].str.contains('-simple-')]
#dat = dat[dat['ResourceComplexity'].str.contains('-monoculture-')]
#dat = dat[dat['TrophicComplexity'].str.contains('-scavenging-')]
#dat = dat[dat['TrophicComplexity'].str.contains('-crossfeeding-')]

dat1 = dat[dat['SpatialComplexity'].str.contains('-chemotaxis-')]
dat2 = dat[dat['SpatialComplexity'].str.contains('-randwalk-')]
dat3 = dat[dat['SpatialComplexity'].str.contains('-none-')]

#-------------------------------------------------------------------------------

#### plot figure ###############################################################
fig = plt.figure()
fs = 8 # fontsize
xlimit = 2000
nlim = 4
lwidth = 1


#### PLOT 1 #################################################################
fig.add_subplot(2, 2, 1)


plt.plot([0, 1], [-10, -10], c = '0.3', ls = '-', lw=2, label='Chemotaxis')
plt.plot([0, 1], [-10, -10], c = 'm', ls = '-', lw=2, label='Run and tumble')
plt.plot([0, 1], [-10, -10], c = 'c', ls = '-', lw=2, label='No active dispersal')

xlab = 'time'
ylab = '% Dormant'

fig = plotlines(fig, dat1, '%Dormant', '0.3', lwidth, nlim)
fig = plotlines(fig, dat2, '%Dormant', 'm', lwidth, nlim)
fig = plotlines(fig, dat3, '%Dormant', 'c', lwidth, nlim)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.xlim(0, xlimit)
#plt.ylim(-2.0, 0.1)
plt.tick_params(axis='both', which='major', labelsize=fs)
plt.legend(bbox_to_anchor=(-0.04, 1.05, 2.48, .2), loc=10, ncol=3, mode="expand",prop={'size':fs+5})


#### PLOT 2 #################################################################
fig.add_subplot(2, 2, 2)

xlab = 'time'
ylab = '$log$'+r'$_{10}$'+'($N$)'

fig = plotlines(fig, dat1, 'N', '0.3', lwidth, nlim)
fig = plotlines(fig, dat2, 'N', 'm', lwidth, nlim)
fig = plotlines(fig, dat3, 'N', 'c', lwidth, nlim)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.xlim(0, xlimit)
#plt.ylim(-0.1, 2.5)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### PLOT 3 #################################################################
fig.add_subplot(2, 2, 3)

xlab = 'time'
ylab = '$log$'+r'$_{10}$'+'($R$)'

fig = plotlines(fig, dat1, 'R', '0.3', lwidth, nlim)
fig = plotlines(fig, dat2, 'R', 'm', lwidth, nlim)
fig = plotlines(fig, dat3, 'R', 'c', lwidth, nlim)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.xlim(0, xlimit)
#plt.ylim(-0.1, 2.5)
plt.tick_params(axis='both', which='major', labelsize=fs)


#### PLOT 4 #################################################################
fig.add_subplot(2, 2, 4)

xlab = 'time'
ylab = '$log$'+r'$_{10}$'+'($Encounters$)'

fig = plotlines(fig, dat1, 'Encounters', '0.3', lwidth, nlim)
fig = plotlines(fig, dat2, 'Encounters', 'm', lwidth, nlim)
fig = plotlines(fig, dat3, 'Encounters', 'c', lwidth, nlim)

plt.ylabel(ylab, fontsize=fs+5)
plt.xlabel(xlab, fontsize=fs+5)
plt.xlim(0, xlimit)
#plt.ylim(0, 1.1)
plt.tick_params(axis='both', which='major', labelsize=fs)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/TimeSeries-Lines.png', dpi=600, bbox_inches = "tight")
plt.close()
