## Model description following the ODD Protocol
The ODD protocol is standard for describing individual-based models (IBMs), see Grimm et al. 2006. ODD stands for Overview, Design concepts, and Details. Here, we descibe our IBM framework according to the ODD protocol.

### Purpose
Our modeling framework simulates life history of individual organisms, the encounter of organisms with resource particles, the emergence of microbial seed banks, and the evolution of traits in spatially explicit environments under stochastic conditions. The goal is to provide a simulation-based framework for examining conditions under which encounter rates have a robust influence on the emergence of seed banks and the growth and activity dynamics of communities. This is done by simulating dimensions of ecological complexity, variation in model parameters and processes, varying physical dynamics, and by allowing realistic life history trade-offs to emerge from random combinations of trait values.

### Entities & their state variables
**Individual organisms** --- Individuals are distinguished by collections of elements within lists. Individuals undergo changes when randomly sampled from lists. Each specific position in the list corresponds to the same individual. For example, in simulating growth, a model chooses an individual from the list of cell quotas (the probability of reproducing is partly determined by endogenous resources). The first position in this list as well in all other lists of individual attributes corresponds to the same individual.

	Example:
	IndIDs = [1,  2, 33, 14]
	Quota =  [0.1, 0.99, 0.14, 0.05]
	Xpos =   [45, 23, 456, 1]
	Ypos =   [765, 87, 21, 34]
	
	The individual with ID of 1 has a cell quota of 0.1 and is located at position x=45, y=765

These are the lists of attributes for individual organisms:

* time each individual spends in the system
* species ID
* individual ID
* 2D spatial location
* endogenous level of a single primary limiting resource (i.e., cell quota)
* species-specific:
	* metabolic maintenance cost
	* maximum dispersal rate
	* maximum growth rate
	* resource use efficiency
	* environmental optima
* individual-level (varies as a normally distributed random variable with the specific-specific value as the mean):
	* metabolic maintenance cost
	* maximum dispersal rate
	* maximum growth rate
	* resource use efficiency
* products of individual-level metabolism


**Species** --- Each species is characterized by the individuals that share a common set of traits, such as maximum growth rate and environmental optima. Species information is stored in Python dictionaries. In this way, if simplex requires the species ID of an individual it will access the species ID list where each element corresponds to a specific individual. But, if simplex requires the maximum growth rate for an individual, then it finds the species ID and then uses that to access the Python dictionary for maximum specific growth rates of species.

	Example:
	IndIDs = [1,  2, 33, 14]
	spIDs = [12, 32, 11, 6]
	MaxGrowthDict = {12: 0.9, 32: 0.5, 11: 0.8, 6: 0.4}
	
	The individual with ID of 1 belongs to species 12, which has a
  theoretical maximum growth rate of 0.9.

The following are the types of species-level information that are stored in Python dictionaries: 

* metabolic maintenance cost
* maximum dispersal rate
* maximum theoretical growth rate
* resource use efficiency for each resource
* environmental optima


**Resources** --- Resources are simulated at two levels (molecules, particles). Resource particles are designated as 'a', 'b', or 'c'. Particles are then joined together to form molecules (e.g., 'aaaa'). Depending on the model, chemical complexity is simulated by introducing bonds (e.g., '-') that must be broken in order to consume resource particles (e.g., 'aa-aa'). The breaking of these bonds is done via random draws of positions along the molecule, where for example an 'aa-aa' is less likely to be split than is an 'a-a-a-a' particle.

Depending on the model, resource molecules can enter the system in clusters. The degree of spatial aggregation of incoming resource molecules is determined by a 2-dimensional Gaussian distribution, where the standard deviation (dispersion) is chosen at random at the start of the model and the mean (x-y position) is chosen at random at each time point.

Each molecule or resource particle is distinguished by information stored in lists, where each position in each list corresponds to the same molecule or independent particle.

	Example:
	resIDs = ['a',  'bb', 'aaaa', 'b']
	Xpos =   [5.6, 3.4, 5.67, 2.0]
	Ypos =   [8.76, 9.8, 3.2, 4.5]
	
	The particle with ID of 0 is 'a' and is 
	located at position x=5.6, y=8.76
	
The following are the types of information stored about each resource (molecule, particle):

* time in the system
* resource type: a, b, c, or d
* size and complexity, e.g., 'a', 'bb', 'cc-cc'
* 2D spatial location

**Physcial barriers** --- These objects are simulated as discrete 2D spatial coordinates that cannot be occupied by any individual entities. The number, size, and location of physical barriers are chosen at random at the start of each model.

###System level state variables
Each run of an IBM begins with random choices for the values of:

* width (5 to 20)
* height (5 to 20)
* basal flow rate (1.0 to 0.0)
* number, size, and location of physical barriers 
* number and direction of environmental gradients
* rate of stochastic fluctuation in basal flow rate
* degree of stochastic fluctuation 
* degree of syncronized flow of individual particles
	* completely in-sync or completely out of sync 

For example:

	from randparams.py:
	
    width = randint(5,100)
    height = randint(5,100)
    
    # number of barriers
    barriers = randint(1,10)

	 # log-series alpha for metacommunity structure
    alpha = np.random.uniform(0.99, 0.999)
    
    reproduction = choice(['fission', 'sexual'])
    speciation = choice(['yes', 'no'])
    
    # size of starting community
    seedCom = choice([10, 100, 1000])
         
    # m = probability of immigration
    m = choice([0.0, 0.0001, 0.0005, 0.001, 0.005])
    ...	

###Spatial and temporal scales
The two general aspects of scale are grain (a.k.a. resolution) and extent (e.g. total area).

**Spatial extent** --- The environment of the IBMs is two dimensional and can vary along each axis from 5 to 20 discrete units. This makes for a potential total extent of 25 to 10,000 discrete patches, each with a grain of 1 square unit.

Note that all individuals and/or resources move in decimal units the limit of which is determined by Python's decimal precision. This means that individual particles can occupy practically infinite locations within patches and likewise, squeeze through barriers (which only occupy integer x-y coordinates).

**Temporal extent** --- Each IBM was run to a state of mean reversion in total abundance, that is, a point where the number of individuals in the system fluctuates around an average value. Greater detail is given under 'Process overview and scheduling'.

**Grain** --- Grain is the smallest unit over which change can happen. For example, as per the original ODD documentation: “One time step represents one year and simulations were run for 100 years.”  In our IBMs, one time step represented a single generation.

### Process overview and scheduling
**Assembly** --- The user runs a core program (i.e., model.py) that chooses random values for system-level state variables including:

* Number of incoming resource molecules:
	*  0 to 90 per time step
* Maximum specific growth rate:
	*  10 - 50% of cell quota
* Maximum maintenance cost for individuals:
	*  1 - 5% of cell quota
* Maximum specific dispersal rate:
	*  1 to 10% of the environment
* Maximum specific resuscitation rate:
	*  probability of 0.1 - 1.0
* Maximum maintenance energy reduction factor:
	*  Going dormant can reduce maintenance energy between 20 and 100 times
* Incoming resource aggregation:
	*  Determined by the standard deviation of a Normal distribution (0.1 - 0.4)

The core program also decides at random whether fluid dynamics will occur and at what rate of flow.

**Core simulation process** --- The IBM begins to run immediately after assembly. After each iteration of resource inflow where varied number of molecules and particles enter the system, each individual is given the chance to consume, grow, reproduce, starve, die, and to disperse towards resources and environmental optima.

**Duration: A run to mean reversion** --- Once assembled, each IBM simulates ecological processes (birth, death, dispersal, growth, consumption, etc.) until the system reaches a point of mean reversion, i.e., the tendency of a system to reverse a directional change. Mean reversion captures whether a system is fluctuating around a mean value. At each time point past 200 generations, the modeling program conducts an Augmented Dickey-Fuller (ADF) test on the previous 100 generations to determine whether a state of mean reversion has been reached. The ADF is well-explained here: https://www.quantstart.com/articles/Basics-of-Statistical-Mean-Reversion-Testing. Once mean reversion is reached (p < 0.05), the model runs for an additional 100 generations. A model is stopped once mean reversion is determined to have occurred. Unless the number of desired models has been reached, the IBM framework simply constructs another model and runs it to mean reversion.


**Simulated life history** --- Our IBMs simulate growth, reproduction, and death via weighted random sampling. This simulate the partly probabilistic and partly deterministic nature of environmental filtering and individual-level interactions. 

*Immigration*
Individuals enter from any point in the environment. Species identities of inflowing propagules are chosen at random from a uniform distribution. Along with a species ID and species-specific maximum rates of resource uptake, active dispersal, resource efficiencies, cell maintenance, an environmental optima, each propagule is given a unique ID and a starting cell quota that represents the state of internal resources. The cell quota also determines the chance of reproduction. 

*Dispersal*
Our IBMs allow individuals to move towards resource particles. This happens by choosing a random sample of resource particles from the environment, and then moving the individual towards the nearest consumable resource particle at a rate determined by its cell quota and maximum dispersal rate.

*Consumption & growth*
Sampled individuals consume resources according to their species-specific maximum rates of uptake and grow according to species-specific resource efficiencies and maintenance costs. Uptake increases the individual cell quotas and decreases ambient resources. Individual cell quotas are then decreased according to a specific physiological maintenance cost. Likewise, faster rates of uptake incurred a larger energetic cost and hence, lower growth efficiency. 

*Reproduction*
Reproduction is limited to asexual reproduction with the possibility of mutation. Individuals reproduce with a probability determined by the combination of mean cell quota and the proportional match to the environmental optima. The cell quota of each resource is evenly divided between the individual and its daughter. The daughter is given a unique individual ID and the species ID of its mother, but is allowed small mutations in individual-level life history parameters.

*Death* 
Individuals sampled at random will die if their cell quota is equal to or less than maintenance energy.

*Emigration* 
Emigration can result from fluid flow, from random spatial movement, or from active dispersal. Individuals and resource particles are considered to have left the local environment when they pass beyond the 2-dimensional edges.

### Initialization and input data

**Initialization**
The model initiates with a random set of values for state-variables, 100 individuals whose species IDs are drawn at random from a uniform distribution. These values are saved, so that a simplex model could be programmed to replicate an analysis. 

**Inputs**
The modeling platform requires no input data other than those parameters that are chosen at random.

### Design concepts

**Ecological complexity** --- Our IBM platform assembles models from random combinations of variables and processes to generate output data that allow the user to test the general influence of particular variables, processes, or dynamics.

*Spatial complexity*
Our IBMs simulate three levels of spatial complexity.

* White noise: Locations of individual organisms and resource molecules changed at each time step in an uncorrelated way. Hence, every organism and resource particle had the same chance of moving to any location within the environment at each time step in the model. This white-noise model created a well-mixed environment with no dispersal limitation. 
* Aggregation and random walks: Resource molecules enter the environment in clusters. Individuals undergo random walks. The average length of the random walk was a species-specific parameter, and larger dispersal distances carried greater energetic costs. The degree of spatial dispersion among incoming resource clusters was chosen at random. 
* Aggregation and chemotaxis: Resource molecules enter in clusters but individuals were capable of sensing resource molecules based on resource density and distance. 

*Trophic complexity*
Our IBMs simulate four levels of trophic complexity, with the last being a combination of two others.

* The first level had only one trophic interaction, which we refer to as a simple "consumer-resource" model. At this level of trophic complexity, all individuals were solely consumers of inflowing resources. 
* The second level of trophic complexity allowed for the consumption of resources contained in dead bacteria (e.g., Rozen et al., 2009), which is a trophic interaction that we referred to as "scavenging". 
* The third level of trophic complexity simulated a situation in which one group of consumer species generated a metabolic byproduct that could be taken up by a second group of consumer species, which in turn generated a byproduct that served as a resource for a third group of species. This situation was meant to simulate conditions that are characteristic of cross-feeding or syntrophy (Pande et al., 2015). A final level of trophic complexity was characterized by a combination scavenging and cross-feeding, which we saw as more ecologically realistic and complex than scavenging or cross-feeding alone.

*Resource complexity* 
Our IBMs simulate three levels of resource complexity. 

* Monoculture: The simplest condition wherein only one type of resource molecule was supplied. These molecules had no chemical complexity and could be consumed without extracellular enzymatic breakdown (Figure 1). Resource molecules were represented by strings of individual particles (e.g., ‘aaaa’) that could be broken down from each end when encountered by organisms. The number, size, and dispersion of resource molecules were chosen at random (Table 1).

* Polyculture: The second level supplied three different types of resources (e.g., ‘aaaa’, ‘bbbb’, ‘cccc’) that could only be used by a specialist consumer (Figure 1). 

* Lock and key: Simulates the structural complexity of resource molecules, i.e., chemical complexity. Consumers must break down resource molecules at specific locations. For example, a molecule would have a hyphen that simulates a chemical bond separating two groups of resource particles (e.g., ‘aaaa-aaaa’). In order to consume a resource particle, i.e., ‘a’, individuals would need to cleave the ‘aaaa-aaaa’ molecule at the bond and then consume an ‘a’ from one of the two resulting ‘aaaa’ molecules. Because molecules are cleaved at random locations, there is a 1/9 chance of cleaving the ‘aaaa-aaaa’ particle and, likewise, a 1/5 chance of cleaving a ‘bb-bb’ particle, and a 2/8 chance of cleaving a ‘bb-bb-bb’ particle (Figure 1). Therefore, a molecule with a proportionately larger number of bonds requires greater time, and hence, energy to break down.

*Note:* The average total number of individual resource particles (e.g., ‘a’, ‘b’, ‘c’) entering the system at a given time and inflow rate was made to be consistent across levels of resource complexity.

**Nutrient limited growth** --- Individual growth and activity is fueled and limited by resources. Individuals cannot grow in the absence of resources. Without resources, individuals eventually go dormant or die.

**Energetic costs** --- All life history processes impose an energetic cost that is manifested as a reduction in cell quota. In some cases these energetic costs are the results of multiplicative outcomes. For example, the energetic cost of actively dispering x units is the product of x and the energetic cost of moving a single unit.

**Emergence** --- Community- and ecosystem-level patterns emerge from the life history and competitive dynamics at the individual level. These dynamics likewise emerge from initially random combinations of species-specific vital rates (e.g. of growth, dispersal, resuscitation, etc.).

**Fluid dynamics** --- IBMs are chosen at random to simulate fluid and non-fluid dynamics. Our IBM framework uses the Lattice-Boltzmann Method (LBM) to simulate fluid flow. An LBM discretizes the environment into a lattice and attaches to each position in the lattice a number of particle densities and velocities for each of nine directions of movement possible in a 2D environment (N, S, E, W, NE, NW, SE, SW, current position).

**Ecological theory** --- Our IBM platform implicitly draws from multiple bodies of ecological theory.

*Life history theory*   
Life history theory seeks to understand aspects of ontogeny, reproduction, life span, etc. as the result of natural selection. Trade-offs play a central role in life history theory, where investing in one trait (e.g., rapid growth) leads to a decrease in another (e.g., growth efficiency). There are several trade-offs that emerge within the models.

* Mobility vs. metabolic maintenance: Greater active dispersal allows organisms to find resources
* Growth rate vs. metabolic maintenance
* Generalist vs. specialist
* R vs. K selection

*Resource-limitation theory*  
Proposes that growth-limiting resources are the primary determinants of competition and community structure. Our IBMs enforce resource limitation by making all biological processes dependent on levels of endogenous resources. Resources are the primary sources of inter- and intra-specific competition among individuals.

*Ecological niche theory*  
Proposes that each ecological species has a unique range envirommental characteristics within which it can experience positive growth. The more similar the niche of two species, the more intense the competition is likely to be between them. While definitions of ecological niches and how to test niche theory has long been argued, the importance of the ecological niche is a central principle of ecological theory. Our IBMs inherently give each species an ecological niche by assigning each species an environmental optima, a specific resource to use which a subset of other species can also grow on, and unique combination of traits that produce life history trade-offs and differential success in specific environments.

*Ecological neutral theory and idiosyncrasy theory*
Our IBMs operate via random sampling and can vary from being completely neutral (all individuals having equal vital rates) to completely idiosyncratic (all individual and species are as different as possible). However, there is little chance that any given model will be 100% neutral or idiosyncratic. One aspect of neutral theory that our modeling adopts is the importance of stochastic life history processes. This process-based stochasticity is simulated via unweighted random draws among individuals. Once drawn from the community at random, individuals undergo life history processes according to probabilities that are effectively weighted by their cell quota and specific-specific vital rates.

**Hypotheses** --- We used the inherent ecological complexity of our modeling framework to test for robust ecological relationships related to our primary hypothesis that microscale properties of resource encounter rates have a stronger influence on growth,  abundance, and activity than macroscale properties of resource inflow and concentration.

**Adaptation** --- Individuals can move towards their environmental optima. Populations can become aggregated in areas that provide 
favorable intersections of species optima. Species canevolve by the action of the environmental filter on subpopulation variation in state variables.

**Objectives** --- Individuals can seek conditions that match them to their environmental optima and can also seek to acquire resources through active searching, i.e., chemotaxis.

**Learning** --- There is no aspect of individual-based learning in our modeling framework.

**Prediction** --- Individuals do not have the ability to anticipate conditions.

**Sensing** --- In models that include chemotaxis, individuals can sense the proximity of usable resources. In models that include fluid dynamics, individuals can move against the direction of flow.

**Interaction** --- At the moment, individuals interact through excluding each other from resources (e.g. preemption) and by producing resources (i.e., cross-feeding). There is currently no communication.

**Stochasticity** --- The occurrence of life history processes are conducted via random sampling. In this way, population and community dynamics result, in part, from demographic stochasticity. Likewise, the emergence of suites of life history traits proceeds from initially random combinations of traits.

**Collectives** --- Individuals belong to species. Species belong to communities. Depending on the model, communities belong to trophic levels.

**Observation** --- Many models should be run to examine trends in the variation generated. The following is recorded for each model:

* Values of randomly chosen state variables
* Mean and variance for:
	* Total abundance, $N$
	* Species richness, $S$
	* Growth rate (per species & per capita) 
	* Maintenance cost (per species & per capita) 
	* Resource efficiency (per species & per capita) 
	* Active dispersal rate (per species & per capita) 
	* Compositional turnover
		* Bray-Curtis
		* Sorensen's
	* Species turnover
		* Whittaker's $\beta$ 
	* Species evenness
		* Smith and Wilson's evenness, $E_{var}$
		* Simpson's evenness, $E_{1/D}$ 
	* Species diversity
		* Shannon's diversity, $H'$ 
		* Simpson's diversity, $D_{1/D}$
	* Dominance
		* Absolute, $N_{max}$
		* Relative, $N_{max}/N$
	* Productivity
		* In number of individuals
		* In amount of cell quota
	* Individual residence time, i.e., amount of time individuals spend in the system

These data are stored in file as csv files and can are directly imported into the Python environment via the scripts provided in the project repository: https://github.com/LennonLab/Micro-Encounter/tree/master/fig-scripts


### Submodels & Equations

**Cell quota model of Droop**
In simplex models, individuals grow according to their amounts of endogenous resources (cell quota).

Droop (1968, 1983) gave a relationship between specific growth rate ($\mu$) and cell quota ($Q$):

$$\mu = \mu'_{m} (1 - k_{q}/Q)$$

where $k_{q}$ is the minimum cell quota needed for life, also referred to as the subsistence quota. $\mu'$ is the 

**Maintenance cost of Pirt**
Pirt (1965) states "The variation, with growth rate, of the yield of organism from the substrate used as energy source is attributed to consumption of energy at a constant rate for cell maintenance." He derives a relationships between the growth yield (biomass), the growth rate, and metabolic maintenance.

simplex models use Pirt's concept of a constant maintenance requirement. simplex also draws from Pirt's simple relation for substrate use:  

$$use(total) = use(maintenance) + use(growth)$$

Respiration and activity without growth is not accounted for.

## Notes on simplex source code
Our models operate primarily on lists in a programmatic way, e.g., quickly sorting lists, and removing and returning an element from lists with very little overhead. Likewise, our models generate and hold a lot of information about all the particles and elements in the system, which can become a computationally intensive task. To this end, the modeling coded is written in Python, an easy to read high-level programming language that has many scientific, plotting, and animation libraries. Python gives greater control over the operating system than data analysis languages (e.g. R, Matlab) that can be comparatively slow at purely computational tasks and can greatly limit the amount of memory held in any data object, and even fail to import large amounts of data. Python can also obtain C-like speeds when implementing certain software, e.g., Cython. The output is a broad array of information held in a single csv file.

## Speed & Memory

Our models do not complete until the time series of total abundance values reaches a state of mean reversion (i.e., stationarity). Because simplex models can range from quickly flowing systems with high disturbance to barely flowing and highly stable but depleted systems, simulations can potentially take several minutes or more to complete. Likewise, the ability to simulate many complex scenarios also allows for very large total abundances, the values of which cannot be predicted *a priori* and can even potentially outstrip a computer's memory. 

I ran simplex on a Mid 2010 MacBook Pro (OS X 10.9.5) with a 2.4 GHz Intel Core 2 Duo processor and 4GB of Memory. This system probably represents a below average capacity for modern personal computers, which for this study, was desirable as simplex should be able to be ran on both personal computers and high capacity remote servers. I present results for time to completion and required memory in the Results.

## Animations

Simplex can generate animations of its models and store them in various image file formats. It does this using the matplotlib animation libary. At the moment, choosing whether to animate or not animate a model is done by commenting out lines of code in model.py. However, this feature will be developed more highly for future versions.

#Results

## Unit tests

simplex passed all units tests for 15 diversity indices, ensuring that each index returns either the correct calculated value, or 'NaN' if given any values that cannot be used (e.g., negative numbers, string characters, empty lists).

## Speed & Memory

**Results from 100 randomly assembled models.** On average, simplex models required 96.60 +- 32.70 seconds to run, had an average total abundance (*N*) of 31,367 +- 879, and required 161 +- 4.3 megabytes of memory. The longest any model took to complete was 49 minutes. This particularly demanding model required 261MB of memory and also had generated the greatest total abundance (*N* = 68,808). The shortest any model took to complete was 5.8 seconds with a total abundance of 0 individuals. The least amount of memory any model required was 90MB. These and other analyses from the 100 randomly assembled models can be run using the SpeedMemory.Rmd file located in the results/TestResults directory.

## Products
### Output data
simplex generates six files as its output data. They are: 

*SimData.csv*--Formatted as an R data frame, where each row is a run of a randomly assembled model, and each column holds a piece of data about the system that was modeled (e.g., flow rate, total abundance, species richness, species turnover, rate of disturbance, etc.).  

*IndRTD.csv*, *ResRTD.csv*, *TracerRTD.csv*--Each line of these three files contains the amount of time that each individual, resource particle, or inert tracer spent in the system. Each line represents the same run of a randomly assembled model.

*RADs.csv*--Each line is a run of a randomly assembled model and contains the rank-abundance vector at the point when the model was stopped.

*Species.csv*--Each line is a run of a randomly assembled model and contains a vector of species labels corresponding to RADs.csv.

All output data were correctly formatted and placed within results/simulated_data/examples directory.
Each of the six output data files was able to be imported into the RStudio environment using the Exploratory.Rmd Rmarkdown file provided in the "GitHub/simplex/results/analyses/Rmd/" path. The user can use the Exploratory.Rmd file to craft a .Rmd file for their own specific analyses.