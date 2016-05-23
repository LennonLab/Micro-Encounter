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

**MeanResourceConcentration and VarResourceConcentration:** Time-averaged mean and variance for the ratio of the number of resource particles to the area of the environment.

**MeanResourceParticles and VarResourceParticles:** Time-averaged mean and variance for the number of resource particles in the system.

**Speciation:** Per capita probability of a point speciation event, the same as modeled in Hubbell's neutral theory.

**MeanPerCapitaGrowth and VarPerCapitaGrowth:** Time-averaged mean values among individuals.

**MeanPerCapitaMaint and VarPerCapitaMaint:** Time-averaged mean values among individuals.

**MeanPerCapitaActiveDispersal and VarPerCapitaActiveDisperal:** Time-averaged mean values among individuals.

**MeanSpecGrowth and VarSpecGrowth:** Time-averaged mean values among species.

**MeanSpecDisp and VarSpecDisp:** Time-averaged mean values among species.

**MeanSpecMaint and VarSpecMaint:** Time-averaged mean values among species.

**MeanSpecDist and VarSpecDist:** Time-averaged mean values among species.

**MeanDormFreq and VarDormFreq:** Time-averaged mean values among species.

**TrophicComplexityLevel:**

**ResourceComplexityLevel:**

**SpatialComplexityLevel:** 
1. White noise, a.k.a, uncorrelated movement. The x-y coordinates of all particles change in a uniform random way. This effectively removes spatial structure from the model, making it analogous to a thoroughly mixed chemostat or batch reactor.  
2. Aggregated resources and random walking individuals. Resources enter in clumps, which vary in size according to the standard deviation of a 2-D Gaussian distribution (see IncomingResAgg below). Individuals undergo correlated random walks, where the positions at subsequent time steps are influenced by the position at the previous time step.  
3. Aggregated resources and individuals that que in on resource particles (chemotaxis of sorts). Individuals sense which consumable particles are nearest and move towards them.


**MeanEncounter and VarEncounter:** Mean and variance in a time-series of encounter values. Here, an encounter is any event where an individual touches a consumable resource particle.

**IncomingResAgg:** The standard deviation used to parameterize a 2-dimensional Gaussian distribution that determines the size of clumps that resource particles enter as.

**MeanIndAgg and VarIndAgg, MeanResAgg and VarResAgg:** Time-averaged mean and variance for aggregation of individual organisms and resource particle. Quantified using Morisita's index of dispersion. A random distribution produces values close to 1. A uniform distribution produces values of 0 or less. Aggregated values produce values between 1 and the number of quadrats/windows/sections/etc.

**RunTime:** Number of generations a model ran. In starting with a small seed (n = 10), models hit mean reversion very quickly.