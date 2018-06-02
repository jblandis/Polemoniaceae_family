# Polemoniaceae_family
Analysis scripts for Polemoniaceae family phylogenetic study

The scripts in the repository allow for analyses carried out for the Polemoniaceae family phylogenetic analyses. Data file for each analysis was a CSV with the first column being species names and the second (and in some cases the third as well) being the data of interest. Only exception is for the HiSSE analyses where the first two columns are the species names (identical) and the thir column the trait of interest. Specifics for each script:

- ASR_Diversitree_Pollinators.R allows for ancestral state reconstruction for a binary trait incorporating differences in diversification
- ASR_phytools.R allows for ancrstral state reconstructions using phytools for both binary and multistate data
- BiSSE_pollinators.R allwows for BisSE analysis with a binary trait, performing both maximum-likelihood test and mcmc analysis for plots
- FiSSE_pollinators.R allows for FiSSE test with a binary trait
- HiSSE_pollinators_deadend.R specifically tests whether selfing is an evolutionary deadend in a HiSSE framework by not allowing transitions from selfing to outcrossing
- HiSSE_pollinators.R tests 25 models in a HiSSE framework to determine which is the preferred model; as well as plotting the preferred model
- Indpendent_Contrast/PGLS.R allows for multiple tests including indendpent contrast and PGLS for corolla length and corolla width; as well as calculating Pagel's lambda and Blomberg's K for each trait
- MuSSE.R allows for diversification analysis of flower color with four possible states, including maximum-likelihood analysis and mcmc analysis, as well as figure plotting
- SIMMAP_Phytools.R performs stochastic mapping for both a binary trait (autogamy versus outcrossing) and multistate (4 possible flower colors) over one tree, a distribution of trees, and when incorporating the transition matrix from the MuSSE analysis.
-Continuous_ASR.R perform ancestral state reconstructions for corolla length and width at opening of corolla tube
-Phylogenetic_PCA_and_clustering.R allows for phylogenetic PCA of length and width, as well as clustering algorithm to see if distinct clades are formed based on corolla length and width.
