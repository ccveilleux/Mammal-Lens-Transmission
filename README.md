**Lens UV Filtering as an Adaptation for Spatial Resolution in Mammals**

Data and R code

**Contents:**
- Mammal_UV_Dataset.csv -- comparative dataset of activity pattern, visual acuity, corneal diameter, eye axial length, relative cornea size, and lens UV transmission metrics for 38 mammalian species. References for data are provided in Supplemental Table S1 in the manuscript. Lens UV transmission data were collected from Douglas and Jeffery (2014).

- mammal_trees_1000.nex -- subset of 1000 phylogenetic trees pruned to the 38 focal species that were used in all phylogenetic regression analyses.

- LensTransmission_Analysis.R -- Full analysis script: data preparation, tree preparation, phylogenetic regressions (Tables 1-3), collinearity diagnostics, and figure generation.

**Requirements:**
All analyses were run in R 4.5.1. Required packages: ape (v.5.8-1), car (v.3.1-5), dplyr (v.1.2.0), geiger (v.2.0.11), ggplot2 (v.4.0.2), phylolm (v.2.6.5), phytools (v.2.5-2).

**Tree provenance:**
The full set of 10,000 mammalian phylogenetic trees is from VertLife.org:
MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_FBDasZhouEtAl_all10k_v2_nexus.trees

From: Upham NS, Esselstyn JA, Jetz W. 2019. Inferring the mammal tree: Species-level sets of
phylogenies for questions in ecology, evolution, and conservation. PLOS Biology 17, e3000494.
(https://doi.org/10.1371/journal.pbio.3000494)

mammal_trees_1000.nex (included here) is a random 1000 tree subsample drawn from that set using set.seed(247) in R, then pruned to the 38 species in this dataset. Code to reproduce the subsampling is provided in the R file.

**Running the analysis:**
The R script is organized into four sections:
1. Prepare Mammalian Dataset -- loads and transforms the csv data.
2. Prepare Trees -- subsamples and prunes the phylogenetic trees. This can be skipped by loading mammal_trees_1000.nex directly, as shown in the script
3. Perform phylogenetic regressions -- runs the models reported in Tables 1-3. Note: within each analysis block, several candidate formulas and data subsets are listed sequentially. To reproduce a specific table result, choose the corresponding formula/subset pair before performing run_phylogenetic_analysis().
   Analysis 1 (AP.3 / AP.2 / AP.4 formulas) --> Table 1 (activity pattern models)
   Analysis 2 (log_VA formulas) --> Table 2 (visual acuity models: All Taxa, Non-Haplorhines, Non-Diurnal Taxa subsets)
   Analysis 3 (Rel.Cornea formulas) --> Table 3 (relative cornea size models, same subset structure as Analysis 2), and Table S2 (log_VA ~ Rel. Cornea validation model)
4. Making Plots -- generates Fig 1 panels and Fig S1.

Expect console warnings during Analysis 2 and 3 indicating that the models are dropping taxa with missing data for the given formula - this is expected behavior because sample sizes differ by predictor and subset, as reported in Tables 2 and 3.

The VIF collinearity check (car::vif on a non-phylogenetic lm) is a diagnostic check only and not one of the reported models in the manuscript tables.

Contact:
Carrie Veilleux (cveill AT midwestern.edu) for questions about the data or code.
