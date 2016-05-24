#Metadata
A file to provide information on simulation data in any of the *SimData.csv files. 

##Columns and descriptions:

**RowID.** A number to corresponding to a particular model. Is not currently used in any analysis. In older data files, this number is not necessarily a unique identifier.

**MeanIndProduction and VarIndProduction:** Time-averaged mean number of individuals produced per time step. If this number is negative, it reflects a net loss.

**ResInflow:** Average number of resource particles entering the system per time step.

**MaxGrowthRate:** Maximum per capita growth rate for any species. This is the highest proportional increase in cell quota that can any individual can attain before maintenance costs are subtracted. A species specific parameter.

**MaxMetMaint:** Maximum value by which the cell quota of any individual is decreased per time step. A species-specific parameter value, the influence of which is independent of the cell quota value. That is, a smaller cell quota does not result in less maintenance.

**MaxActiveDispersal:** Maximum distance an individual can disperse; a species species parameter. The value relates to a portion of the landscape. E.g., able to disperse (at most) 0.1% of the landscape in a give time step. Note: spatial grain and extent of the environment is rather arbitrary.

**LogSeriesA:** The parameter value for the log-series distribution. The log-series is a successful model of (meta)commuity structure. Each model starts with a random sample of individuals whose species ID's are drawn from a log-series distribution.

**StartingSeed:** The number of individuals a model is started with. The larger the number, the longer a model takes to escape the "burn in" period. For fast runs that yeild a few thousand models in one day, a seed of 10 works well. A seed of 100 might produce a few hundred models in one day, but will often result in value of total individual abundance that can bog down a humble laptop.

**width and height:** Width and height of the environment; an arbitrary extent and grain.

**MeanTotalAbundance and VarTotalAbundance:** Time-averaged mean and variance for the number of individuals in the system.

**Immigration rate:** Number of individuals entering the system from outside, per time step.

**MeanResourceConcentration and VarResourceConcentration:** Mean and variance in a time-series of values for the ratio of the number of resource particles to the area of the environment.

**MeanResourceParticles and VarResourceParticles:** Mean and variance in a time-series of values quantifying the number of resource particles in the system.

**Speciation:** Per capita probability of a point speciation event, the same as modeled in Hubbell's neutral theory.

**MeanPerCapitaGrowth and VarPerCapitaGrowth:** Mean and variance in a time-series of values quantifying per capita growth rates.

**MeanPerCapitaMaint and VarPerCapitaMaint:** Mean and variance in a time-series of per capita metabolic maintenance values.

**MeanPerCapitaActiveDispersal and VarPerCapitaActiveDisperal:** Mean and variance in a time-series of per capita active dispersal rates.

**MeanSpecGrowth and VarSpecGrowth:** Mean and variance in a time-series of specific growth rates.

**MeanSpecDisp and VarSpecDisp:** Mean and variance in a time-series of dispersal rate values.

**MeanSpecMaint and VarSpecMaint:** Mean and variance in a time-series of species specific metabolic maintenance values.

**MeanSpecDist and VarSpecDist:** Mean and variance in a time-series of values quantifying the species specific average distance to resource particles.

**MeanDormFreq and VarDormFreq:** Mean and variance in a time-series of values quantifying the fraction of total abundance that is dormant.

**TrophicComplexityLevel:** One of three levels of trophic complexity.  
1. A single trophic level. Species do not produce resource for other species. There is no scavenging or recycling of dead biomass.  
2. Three trophic levels. Resource B is produced as a by-product of metabolizing resource A. Resource C, in turn, is produced as a by-product of metabolizing resource B.  
3. Recycling of dead individuals. This level only differs from level #1 in the recycling. Dead individuals become available to all species as resource "D".

**ResourceComplexityLevel:** One of three levels of resource complexity.  
1. All resource particles are a's.  
2. Resource particles are a's, b's, or c's.  
3. Resource particles are a's, b's, or c's, and all particles have a structure complexity.

**SpatialComplexityLevel:** One of three levels of spatial complexity.  
1. White noise, a.k.a, uncorrelated movement. The x-y coordinates of all particles change in a uniform random way. This effectively removes spatial structure from the model, making it analogous to a thoroughly mixed chemostat or batch reactor.  
2. Aggregated resources and random walking individuals. Resources enter in clumps, which vary in size according to the standard deviation of a 2-D Gaussian distribution (see IncomingResAgg below). Individuals undergo correlated random walks, where the positions at subsequent time steps are influenced by the position at the previous time step.  
3. Aggregated resources and individuals that que in on resource particles (chemotaxis of sorts). Individuals sense which consumable particles are nearest and move towards them.

**MeanEncounter and VarEncounter:** Mean and variance in a time-series of encounter values. Here, an encounter is any event where an individual touches a consumable resource particle.

**IncomingResAgg:** The standard deviation used to parameterize a 2-dimensional Gaussian distribution that determines the size of clumps that resource particles enter as.

**MeanIndAgg and VarIndAgg, MeanResAgg and VarResAgg:** Mean and variance in a time-series of aggregation values pertaining to individual organisms and resource particles. Quantified using Morisita's index of dispersion. A random distribution produces values close to 1. A uniform distribution produces values of 0 or less. Aggregated values produce values between 1 and the number of quadrats/windows/sections/etc.

**RunTime:** Number of generations a model ran. In starting with a small seed (n = 10), models hit mean reversion very quickly.