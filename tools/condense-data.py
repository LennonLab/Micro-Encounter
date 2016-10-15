from __future__ import division
import linecache
from math import isnan
import numpy as np
from numpy import mean
import pandas as pd
import sys
import os

mydir = os.path.expanduser("~/GitHub/Micro-Encounter")

sys.path.append(mydir + "/model/col_labels")
labels = linecache.getline(mydir + '/model/col_labels/condensed-data.txt', 1)
with open(mydir + '/results/simulated_data/SimData.csv', 'w+') as text_file:
    text_file.write(labels)


df = pd.read_csv(mydir + '/results/simulated_data/2016_09_18_SimData.csv')
#df = df.convert_objects(convert_numeric=True).dropna()
#df = df[df['sim'] > 500]

#-------------------------DATA TRANSFORMATIONS -----------------------
df2 = pd.DataFrame({'Encounters' : np.log10(df['Encounters'].groupby(df['sim']).mean())})
df2 = df2[np.isfinite(df2['Encounters'])]

df2['SpatialComplexity'] = df['SpatialComplexity'].groupby(df['sim']).unique()
df2['SpatialComplexity'] = [df2['SpatialComplexity'][i][0] for i in df2['SpatialComplexity'].keys()]

df2['TrophicComplexity'] = df['TrophicComplexity'].groupby(df['sim']).unique()
df2['TrophicComplexity'] = [df2['TrophicComplexity'][i][0] for i in df2['TrophicComplexity'].keys()]

df2['ResourceComplexity'] = df['ResourceComplexity'].groupby(df['sim']).unique()
df2['ResourceComplexity'] = [df2['ResourceComplexity'][i][0] for i in df2['ResourceComplexity'].keys()]

df2['NumberDead'] = df['numDead'].groupby(df['sim']).mean()
df2['TotalAbundance'] = df['N'].groupby(df['sim']).mean()
df2['DormantN'] = df['DormantN'].groupby(df['sim']).mean()
df2['%Dormant'] = df2['DormantN']/df2['TotalAbundance']
df2['Productivity'] = df['PRODI'].groupby(df['sim']).mean()
df2['ActiveN'] = df2['TotalAbundance'] - df2['DormantN']
df2['TotalResources'] = df['R'].groupby(df['sim']).mean()
df2['ResourceInflow'] = df['ResInflow'].groupby(df['sim']).mean()

# TRAITS
df2['MeanCellQuota'] = df['MeanCellQuota'].groupby(df['sim']).mean()
df2['MaxGrowth'] = df['MaxGrowth'].groupby(df['sim']).mean()
df2['MaxMaint'] = df['MaxMaint'].groupby(df['sim']).mean()
df2['MaxDispersal'] = df['MaxDispersal'].groupby(df['sim']).mean()
df2['MaxRPF'] = df['MaxRPF'].groupby(df['sim']).mean()
df2['MaxMaintFactor'] = df['MaxMainFactor'].groupby(df['sim']).mean()
df2['SpeciesSpecificDispersal'] = df['SpeciesDisp'].groupby(df['sim']).mean()
df2['SpeciesSpecificMaintenance'] = df['SpeciesMaint'].groupby(df['sim']).mean()
df2['SpeciesSpecificGrowth'] = df['SpecificGrowth'].groupby(df['sim']).mean()
df2['PerCapitaGrowth'] = df['PerCapitaGrowth'].groupby(df['sim']).mean()
df2['PerCapitaMaint'] = df['PerCapitaMaint'].groupby(df['sim']).mean()
df2['PerCapitaDispersal'] = df['PerCapitaDisp'].groupby(df['sim']).mean()

path = mydir + '/results/simulated_data/SimData_condensed.csv'
df2.to_csv(path, sep=',')
