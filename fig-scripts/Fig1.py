from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys

#import statsmodels.stats.api as sms
#import statsmodels.api as sm
import statsmodels.formula.api as smf
#from statsmodels.sandbox.regression.predstd import wls_prediction_std
from statsmodels.stats.outliers_influence import summary_table


mydir = os.path.expanduser('~/Desktop/MicroEncounter')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")

#dat = pd.read_csv(mydir + '/results/simulated_data/examples/2015_01_19/SimData.csv')
dat = pd.read_csv(mydir + '/SimData/SimData.csv')

#print 'here'
#sys.exit()

dat = dat[np.isfinite(dat['resource.concentration'])]
dat = dat[np.isfinite(dat['resource.particles'])]
dat = dat[np.isfinite(dat['dorm.freq'])]
dat = dat[np.isfinite(dat['avg.dist'])]

dat = dat[dat['avg.dist'] > 0]
#dat = dat[dat['resource.concentration'] > 0]
#dat = dat[dat['resource.particles'] > 0]
#dat = dat[dat['barriers'] == 0]

dat['dorm'] = dat['dorm.freq']
dat['res_conc'] = np.log10(dat['resource.concentration'])
dat['res_N'] = np.log10(dat['resource.particles'])
dat['avg_dist'] = dat['avg.dist']

Dorm = dat['dorm'].tolist()
ResConc = dat['res_conc'].tolist()
NRes = dat['res_N'].tolist()
AvgDist = dat['avg_dist'].tolist()

#### plot figure ###############################################################
fs = 8 # fontsize
fig = plt.figure()

#### Itau vs. Tau #################################################################
fig.add_subplot(2, 2, 1)

xlab = 'Percent Dormant'
#f2 = smf.ols('res_conc ~ dorm + I(dorm ** 2.0)', dat).fit()
f2 = smf.ols('res_conc ~ dorm', dat).fit()
#print f2.summary(),'\n\n'
#f2 = smf.ols('Itau ~ tau', d).fit()

a, b  =  f2.params
p1, p2 = f2.pvalues
#a, b, c =  f2.params
#p1, p2, p3 = f2.pvalues
r2 = round(f2.rsquared, 2)

st, data, ss2 = summary_table(f2, alpha=0.05)
fitted = data[:,2]
pred_mean_se = data[:,3]
pred_mean_ci_low, pred_mean_ci_upp = data[:,4:6].T
pred_ci_low, pred_ci_upp = data[:,6:8].T

dorm2, fitted, pred_ci_low, pred_ci_upp, pred_mean_ci_low, pred_mean_ci_upp = zip(*sorted(zip(Dorm, fitted, pred_ci_low, pred_ci_upp, pred_mean_ci_low, pred_mean_ci_upp)))
plt.scatter(dorm2, ResConc, color = 'b', alpha = 0.4 , s = 5, linewidths = 0.0, edgecolor = 'k')
plt.fill_between(dorm2, pred_ci_low, pred_ci_upp, color='r', lw=0.0, alpha=0.1)
plt.fill_between(dorm2, pred_mean_ci_low, pred_mean_ci_upp, color='r', lw=0.0, alpha=0.3)
plt.plot(dorm2, fitted, color = 'r', ls='--', lw=0.5, alpha=0.9)

ylab = 'Resource\nconcentration'
plt.ylabel(ylab, fontsize=fs+5)
plt.xlim(0, 1)
#plt.ylim(0, 7)
plt.xlabel(xlab, fontsize=fs+5)
plt.tick_params(axis='both', which='major', labelsize=fs)

#plt.text(4, 0.5,  r'$r^2$' + '=' +str(r2)+', '+r'$p$' + '=' +str(round(p2, 3)), fontsize=fs+1, color='k')



#### Number of resource particles vs. % Dormant ################################

fig.add_subplot(2, 2, 2)
#f2 = smf.ols('res_N ~ dorm + I(dorm ** 2.0)', dat).fit()
f2 = smf.ols('res_N ~ dorm', dat).fit()
#print f2.summary(),'\n\n'
#f2 = smf.ols('Ptau ~ tau', d).fit()
#print f2.summary()

a, b  =  f2.params
p1, p2 = f2.pvalues
#a, b, c =  f2.params
#p1, p2, p3 = f2.pvalues
r2 = round(f2.rsquared, 2)

st, data, ss2 = summary_table(f2, alpha=0.05)
fitted = data[:,2]
pred_mean_se = data[:,3]
pred_mean_ci_low, pred_mean_ci_upp = data[:,4:6].T
pred_ci_low, pred_ci_upp = data[:,6:8].T

dorm2, fitted, pred_ci_low, pred_ci_upp, pred_mean_ci_low, pred_mean_ci_upp = zip(*sorted(zip(Dorm, fitted, pred_ci_low, pred_ci_upp, pred_mean_ci_low, pred_mean_ci_upp)))
plt.scatter(dorm2, NRes, color = 'b', alpha = 0.4 , s = 5, linewidths = 0.0, edgecolor = 'k')
plt.fill_between(dorm2, pred_ci_low, pred_ci_upp, color='r', lw=0.0, alpha=0.1)
plt.fill_between(dorm2, pred_mean_ci_low, pred_mean_ci_upp, color='r', lw=0.0, alpha=0.3)
plt.plot(dorm2, fitted, color = 'r', ls='--', lw=0.5, alpha=0.9)


plt.ylabel('Resource\nparticles', fontsize=fs+5)
plt.xlim(0, 1)
#plt.ylim(0, 3)
plt.xlabel(xlab, fontsize=fs+5)
plt.tick_params(axis='both', which='major', labelsize=fs)

#plt.text(3.5, 0.5,  r'$r^2$' + '=' +str(r2)+', '+r'$p$' + '=' +str(round(p2, 3)), fontsize=fs+1, color='k')

####  AvgDist vs. % Dormant #################################################################
fig.add_subplot(2, 2, 3)

#f2 = smf.ols('avg_dist ~ dorm + I(dorm ** 2.0)', dat).fit()
f2 = smf.ols('avg_dist ~ dorm', dat).fit()
print f2.summary(),'\n'
#sys.exit()

a, b  =  f2.params
p1, p2 = f2.pvalues
#a, b, c =  f2.params
#p1, p2, p3 = f2.pvalues
r2 = round(f2.rsquared, 2)


st, data, ss2 = summary_table(f2, alpha=0.05)
fitted = data[:,2]
pred_mean_se = data[:,3]
pred_mean_ci_low, pred_mean_ci_upp = data[:,4:6].T
pred_ci_low, pred_ci_upp = data[:,6:8].T

dorm2, fitted, pred_ci_low, pred_ci_upp, pred_mean_ci_low, pred_mean_ci_upp = zip(*sorted(zip(Dorm, fitted, pred_ci_low, pred_ci_upp, pred_mean_ci_low, pred_mean_ci_upp)))
plt.scatter(dorm2, AvgDist, color = 'b', alpha = 0.4 , s = 5, linewidths = 0.0, edgecolor = 'k')
plt.fill_between(dorm2, pred_ci_low, pred_ci_upp, color='r', lw=0.0, alpha=0.1)
plt.fill_between(dorm2, pred_mean_ci_low, pred_mean_ci_upp, color='r', lw=0.0, alpha=0.3)
plt.plot(dorm2, fitted, color = 'r', ls='--', lw=0.5, alpha=0.9)

ylab = 'Mean distance'
plt.ylabel(ylab, fontsize=fs+5)
plt.xlim(0, 1)
plt.ylim(0, 10)
plt.xlabel(xlab, fontsize=fs+3)
plt.tick_params(axis='both', which='major', labelsize=fs)

#plt.text(0.1, 2.2,  r'$r^2$' + '=' +str(r2)+', '+r'$p$' + '=' +str(round(p2,3)), fontsize=fs+1, color='k')


#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir2 + 'Desktop/MicroEncounter/Figures/Fig1.png', dpi=600, bbox_inches = "tight")
#plt.show()
