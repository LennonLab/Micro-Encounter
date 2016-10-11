from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from random import randint, sample
import os
import sys
import statsmodels.api as sm


def randcolor():
    c1 = randint(0,255)
    c2 = randint(0,255)
    c3 = randint(0,255)

    clr = '#%02x%02x%02x' % (c1, c2, c3)
    return clr


def plotlines(fig, dat, datlabel, clr, lwidth, sims):

    for sim in sims:
        d = dat[dat['sim'] == sim]
        x = d['ct'].tolist()
        y = np.log10(d[datlabel]).tolist()

        lowess = sm.nonparametric.lowess(y, x, frac=0.01)
        plt.plot(lowess[:, 0], lowess[:, 1], c = clr, ls = '-', lw=lwidth)

    return fig

mydir = os.path.expanduser('~/GitHub/Micro-Encounter')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")
df = pd.read_csv(mydir + '/results/simulated_data/SimData-TSA.csv')
#df = df.convert_objects(convert_numeric=True).dropna()

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
nlim = 1
lwidth = 1


#### PLOT 1 #################################################################
fig.add_subplot(2, 2, 1)

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


#### PLOT 2 #################################################################
fig.add_subplot(2, 2, 2)

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
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/TimeSeriesLines-BoomBust.png', dpi=600, bbox_inches = "tight")
plt.close()
