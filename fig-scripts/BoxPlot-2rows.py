from __future__ import division
import  matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
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

    if clr == 'k': colors = ['0.6']*len(bp['boxes'])
    elif clr == 'b': colors = ['c']*len(bp['boxes'])
    elif clr == 'purple': colors = ['m']*len(bp['boxes'])

    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)


#df = pd.read_csv(mydir + '/results/simulated_data/SimData.csv')
df = pd.read_csv(mydir + '/results/simulated_data/2016_09_18_SimData.csv')
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

#df2 = df2[~df2['SpatialComplexity'].str.contains('-wellmixed-')] # <---------------

#-------------------------END DATA FILTERS & TRANSFORMATIONS -------------------

#### plot figure 1 #############################################################

def boxplotrow(i, df2, fig):

    ax = fig.add_subplot(2, 2, i)

    if i == 1:
        df2 = df2[df2['SpatialComplexity'].str.contains('-wellmixed-')]
    else:
        df2 = df2[df2['SpatialComplexity'].str.contains('-brownian-')]

    if i == 1:
        p1 = Rectangle((0, 0), 1, 1, fc='0.6', ec='k')  # DeepSkyBlue', ec='b'
        p2 = Rectangle((0, 0), 1, 1, fc='m', ec='purple')  # 'y', ec='orange'
        p3 = Rectangle((0, 0), 1, 1, fc='c', ec='b') # 'Pink', ec='r'
        plt.legend([p1, p2, p3], ['Chemotaxis', 'Run and tumble', 'No dispersal'], bbox_to_anchor=(-0.04, 1.05, 2.48, .2), loc=10, ncol=3, mode="expand",prop={'size':fs+5})


    # CHEMOTAXIS
    df1 = df2[np.isfinite(df2['Encounters'])]
    df1 = df1[df1['SpatialComplexity'].str.contains('-chemotaxis-')] # <---------------
    print df1.shape

    dat1 = df1[df1['ResourceComplexity'].str.contains('-monoculture-')]
    dat1 = dat1[~dat1['ResourceComplexity'].str.contains('-lockandkey-')]
    encounters1 = dat1['Encounters']
    print dat1.shape

    dat2 = df1[df1['ResourceComplexity'].str.contains('-polyculture-')]
    dat2 = dat2[~dat2['ResourceComplexity'].str.contains('-lockandkey-')]
    encounters2 = dat2['Encounters']
    print dat2.shape

    dat3 = df1[df1['ResourceComplexity'].str.contains('-polyculture--polymers--lockandkey-')]
    encounters3 = dat3['Encounters']
    print dat3.shape,'\n'

    data_to_plot = [encounters1, encounters2, encounters3]
    bp = ax.boxplot(data_to_plot, showfliers=False, patch_artist=True)

    clr='k'
    setBoxColors(bp, clr)

    # RUN AND TUMBLE

    df1 = df2[np.isfinite(df2['Encounters'])]
    df1 = df1[df1['SpatialComplexity'].str.contains('-randwalk-')] # <---------------

    dat1 = df1[df1['ResourceComplexity'].str.contains('-monoculture--polymers--simple-')]
    encounters1 = dat1['Encounters']
    print dat1.shape

    dat2 = df1[df1['ResourceComplexity'].str.contains('-polyculture--polymers--simple-')]
    encounters2 = dat2['Encounters']
    print dat2.shape

    dat3 = df1[df1['ResourceComplexity'].str.contains('-polyculture--polymers--lockandkey-')]
    encounters3 = dat3['Encounters']
    print dat3.shape,'\n'

    data_to_plot = [encounters1, encounters2, encounters3]
    bp = ax.boxplot(data_to_plot, showfliers=False, patch_artist=True)

    clr='purple'
    setBoxColors(bp, clr)


    df1 = df2[np.isfinite(df2['Encounters'])]
    df1 = df1[df1['SpatialComplexity'].str.contains('-none-')] # <---------------

    dat1 = df1[df1['ResourceComplexity'].str.contains('-monoculture--polymers--simple-')]
    encounters1 = dat1['Encounters']
    print dat1.shape

    dat2 = df1[df1['ResourceComplexity'].str.contains('-polyculture--polymers--simple-')]
    encounters2 = dat2['Encounters']
    print dat2.shape

    dat3 = df1[df1['ResourceComplexity'].str.contains('-polyculture--polymers--lockandkey-')]
    encounters3 = dat3['Encounters']
    print dat3.shape,'\n'

    data_to_plot = [encounters1, encounters2, encounters3]
    bp = ax.boxplot(data_to_plot, showfliers=False, patch_artist=True)

    clr='b'
    setBoxColors(bp, clr)

    ax.set_xticklabels(['Monoculture', 'Polyculture', 'Lock & Key'], rotation=45)
    plt.tick_params(axis='both', which='major', labelsize=fs)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_ylabel('Avg encounter, '+ 'log'+r'$_{10}$', fontsize=fs+2)

    if i == 1:
        plt.text(-0.75, 1.3, 'Well-mixed', fontsize=fs+12, rotation=90)
    else:
        plt.text(-0.75, 1.2, 'Structured', fontsize=fs+12, rotation=90)


    #### plot figure 2 #############################################################

    ax = fig.add_subplot(2, 2, i+1)

    df1 = df2[np.isfinite(df2['Encounters'])]
    df1 = df1[df1['SpatialComplexity'].str.contains('-chemotaxis-')]

    dat1 = df1[df1['TrophicComplexity'].str.contains('-none-')]
    encounters1 = dat1['Encounters']
    print dat1.shape

    dat2 = df1[df1['TrophicComplexity'].str.contains('-crossfeeding-')]
    encounters2 = dat2['Encounters']
    print dat2.shape

    dat3 = df1[df1['TrophicComplexity'].str.contains('-scavenging-')]
    encounters3 = dat3['Encounters']
    print dat3.shape

    dat4 = df1[df1['TrophicComplexity'].str.contains('-crossfeeding--scavenging-')]
    encounters4 = dat4['Encounters']
    print dat4.shape,'\n'

    data_to_plot = [encounters1, encounters2, encounters3, encounters4]
    bp = ax.boxplot(data_to_plot, showfliers=False, patch_artist=True)

    clr='k'
    setBoxColors(bp, clr)

    df1 = df2[np.isfinite(df2['Encounters'])]
    df1 = df1[df1['SpatialComplexity'].str.contains('-randwalk-')]

    dat1 = df1[df1['TrophicComplexity'].str.contains('-none-')]
    encounters1 = dat1['Encounters']
    print dat1.shape

    dat2 = df1[df1['TrophicComplexity'].str.contains('-crossfeeding-')]
    encounters2 = dat2['Encounters']
    print dat2.shape

    dat3 = df1[df1['TrophicComplexity'].str.contains('-scavenging-')]
    encounters3 = dat3['Encounters']
    print dat3.shape

    dat4 = df1[df1['TrophicComplexity'].str.contains('-crossfeeding--scavenging-')]
    encounters4 = dat4['Encounters']
    print dat4.shape,'\n'

    data_to_plot = [encounters1, encounters2, encounters3, encounters4]
    bp = ax.boxplot(data_to_plot, showfliers=False, patch_artist=True)

    clr='purple'
    setBoxColors(bp, clr)


    df1 = df2[np.isfinite(df2['Encounters'])]
    df1 = df1[df1['SpatialComplexity'].str.contains('-none-')] # <---------------

    dat1 = df1[df1['TrophicComplexity'].str.contains('-none-')]
    encounters1 = dat1['Encounters']
    print dat1.shape

    dat2 = df1[df1['TrophicComplexity'].str.contains('-crossfeeding-')]
    encounters2 = dat2['Encounters']
    print dat2.shape

    dat3 = df1[df1['TrophicComplexity'].str.contains('-scavenging-')]
    encounters3 = dat3['Encounters']
    print dat3.shape

    dat4 = df1[df1['TrophicComplexity'].str.contains('-crossfeeding--scavenging-')]
    encounters4 = dat4['Encounters']
    print dat4.shape,'\n'

    data_to_plot = [encounters1, encounters2, encounters3, encounters4]
    bp = ax.boxplot(data_to_plot, showfliers=False, patch_artist=True)

    clr='b'
    setBoxColors(bp, clr)

    ax.set_xticklabels(['Simple\ncommunity', 'Crossfeeding', 'Scavenging', 'Crossfeeding\nScavending'], rotation=45)
    plt.tick_params(axis='both', which='major', labelsize=fs)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_ylabel('Avg encounter, '+ 'log'+r'$_{10}$', fontsize=fs+2)

    return fig


fs = 8
fig = plt.figure()

i = 1
fig = boxplotrow(i, df2, fig)

i = 3
fig = boxplotrow(i, df2, fig)
#### Final Format and Save #####################################################

plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/BoxPlot-2x2.png', dpi=600, bbox_inches = "tight")
#plt.show()
#plt.close()
