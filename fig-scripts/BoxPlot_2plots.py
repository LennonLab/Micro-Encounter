from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys

mydir = os.path.expanduser('~/GitHub/Micro-Encounter')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")


def setBoxColors(bp, clr):
    plt.setp(bp['boxes'], color=clr)
    plt.setp(bp['caps'], color=clr)
    plt.setp(bp['whiskers'], color=clr)
    plt.setp(bp['medians'], color=clr)

    if clr == 'r': colors = ['Pink']*len(bp['boxes'])
    elif clr == 'b': colors = ['DeepSkyBlue']*len(bp['boxes'])

    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)



df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')
#df = df.convert_objects(convert_numeric=True).dropna()

#-------------------------DATA FILTERS & TRANSFORMATIONS -----------------------

#http://pandas.pydata.org/pandas-docs/stable/generated/pandas.Series.html
#http://chrisalbon.com/python/pandas_apply_operations_to_groups.html
#http://pandas.pydata.org/pandas-docs/stable/merging.html
#http://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.groupby.html
#http://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.to_csv.html
#http://stackoverflow.com/questions/13411544/delete-column-from-pandas-dataframe

#df2['Encounters'] = np.log10(df['Encounters'].groupby(df['sim']).mean())
#df2 = pd.DataFrame({'Encounters' : df['Encounters'].groupby(df['sim']).mean()})
df2 = pd.DataFrame({'Encounters' : np.log10(df['Encounters'].groupby(df['sim']).mean())})

df2['SpatialComplexity'] = df['SpatialComplexity'].groupby(df['sim']).unique()
df2['SpatialComplexity'] = [df2['SpatialComplexity'][i][0] for i in df2['SpatialComplexity'].keys()]

df2['TrophicComplexity'] = df['TrophicComplexity'].groupby(df['sim']).unique()
df2['TrophicComplexity'] = [df2['TrophicComplexity'][i][0] for i in df2['TrophicComplexity'].keys()]

df2['ResourceComplexity'] = df['ResourceComplexity'].groupby(df['sim']).unique()
df2['ResourceComplexity'] = [df2['ResourceComplexity'][i][0] for i in df2['ResourceComplexity'].keys()]


#-------------------------END DATA FILTERS & TRANSFORMATIONS -------------------


#### plot figure 1 #############################################################

fs = 8
fig = plt.figure()
ax = fig.add_subplot(2, 2, 1)

df1 = df2[np.isfinite(df2['Encounters'])]
df1 = df1[df1['SpatialComplexity'].str.contains('-chemotaxis-')] # <---------------

dat1 = df1[df1['ResourceComplexity'].str.contains('-monoculture--polymers--simple-')]
encounters1 = dat1['Encounters']

dat2 = df1[df1['ResourceComplexity'].str.contains('-polyculture--polymers--simple-')]
encounters2 = dat2['Encounters']

dat3 = df1[df1['ResourceComplexity'].str.contains('-polyculture--polymers--lockandkey-')]
encounters3 = dat3['Encounters']

data_to_plot = [encounters1, encounters2, encounters3]
bp = ax.boxplot(data_to_plot, showfliers=False, patch_artist=True)

clr='b'
setBoxColors(bp, clr)

df1 = df2[np.isfinite(df2['Encounters'])]
df1 = df1[~df1['SpatialComplexity'].str.contains('-chemotaxis-')] # <---------------

dat1 = df1[df1['ResourceComplexity'].str.contains('-monoculture--polymers--simple-')]
encounters1 = dat1['Encounters']

dat2 = df1[df1['ResourceComplexity'].str.contains('-polyculture--polymers--simple-')]
encounters2 = dat2['Encounters']

dat3 = df1[df1['ResourceComplexity'].str.contains('-polyculture--polymers--lockandkey-')]
encounters3 = dat3['Encounters']

data_to_plot = [encounters1, encounters2, encounters3]
bp = ax.boxplot(data_to_plot, showfliers=False, patch_artist=True)

clr='r'
setBoxColors(bp, clr)

ax.set_xticklabels(['Monoculture', 'Polyculture', 'Lock & Key'], rotation=45)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.set_ylabel('Avg encounter, '+ 'log'+r'$_{10}$', fontsize=fs+2)


#### plot figure 2 #############################################################

ax = fig.add_subplot(2, 2, 2)

df1 = df2[np.isfinite(df2['Encounters'])]
df1 = df1[df1['SpatialComplexity'].str.contains('-chemotaxis-')] # <---------------

dat1 = df1[df1['TrophicComplexity'].str.contains('-none-')]
encounters1 = dat1['Encounters']

dat2 = df1[df1['TrophicComplexity'].str.contains('-crossfeeding-')]
encounters2 = dat2['Encounters']

dat3 = df1[df1['TrophicComplexity'].str.contains('-scavenging-')]
encounters3 = dat3['Encounters']

dat4 = df1[df1['TrophicComplexity'].str.contains('-crossfeeding--scavenging-')]
encounters4 = dat4['Encounters']

data_to_plot = [encounters1, encounters2, encounters3, encounters4]
bp = ax.boxplot(data_to_plot, showfliers=False, patch_artist=True)

clr='b'
setBoxColors(bp, clr)

df1 = df2[np.isfinite(df2['Encounters'])]
df1 = df1[~df1['SpatialComplexity'].str.contains('-chemotaxis-')] # <---------------

dat1 = df1[df1['TrophicComplexity'].str.contains('-none-')]
encounters1 = dat1['Encounters']

dat2 = df1[df1['TrophicComplexity'].str.contains('-crossfeeding-')]
encounters2 = dat2['Encounters']

dat3 = df1[df1['TrophicComplexity'].str.contains('-scavenging-')]
encounters3 = dat3['Encounters']

dat4 = df1[df1['TrophicComplexity'].str.contains('-crossfeeding--scavenging-')]
encounters4 = dat4['Encounters']

data_to_plot = [encounters1, encounters2, encounters3, encounters4]
bp = ax.boxplot(data_to_plot, showfliers=False, patch_artist=True)

clr='r'
setBoxColors(bp, clr)

ax.set_xticklabels(['Simple\ncommunity', 'Crossfeeding', 'Scavenging', 'Crossfeeding\nScavending'], rotation=45)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.set_ylabel('Avg encounter, '+ 'log'+r'$_{10}$', fontsize=fs+2)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/BoxPlot_2plots.png', dpi=600, bbox_inches = "tight")
#plt.show()
#plt.close()
