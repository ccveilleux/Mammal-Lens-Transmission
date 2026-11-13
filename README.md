**Lens UV Filtering as an Adaptation for Spatial Resolution in Mammals**

Data and R code

**Contents:**
- Mammal_UVA_Dataset.csv -- comparative dataset of activity pattern, visual acuity, corneal diameter, eye axial length, relative cornea size, path length (lens thickness), and lens UV transmission metrics for 47 mammalian species. References for data are provided in Supplemental Table S1 in the manuscript. Lens UV transmission data were collected from Douglas and Jeffery (2014), McDowell et al. (2024), and digitized from hamster (Mesocricetus auratus) lens transmission curve in Hut et al. (2000). Note that McDowell et al. (2024) values were only used for species not measured by Douglas and Jeffery (2014) and represent a subset of the species published in McDowell et al.'s Supplementary Data. Species were only included if the original source studies measured lens transmission across the full UVA range (315-400 nm) on intact lenses (i.e., not compressed). Six species were excluded because either the methods were undescribed (Tupaia belangeri, Sminthopsis crassicaudata), they were only measured in the visible spectrum (400+nm) (Callospermophilus lateralis, Ictidomys mexicanus, Sciurus griseus), or the lenses were compressed during measurement (Meriones unguiculatus). Additionally, the seal (Cystophora cristata) was excluded due to having an aquatic-adapted eye, the European mole (Talpa europaea) for a cellular and permanently cloudy lens, and McDowell's mouse (Mus musculus) measurement because I prioritized Douglas and Jeffery's (2014) mouse measurement.
   
Columns: Species, Order, Path.Length (mm), PropUVA (proportion of UVA transmitted by lens), 50.T.nm (λT0.5, wavelength (nm) at which 50% of incident light is transmitted   by the lens), VA (visual acuity, cpd), Acuity.Type (A-anatomical, B-behavioral), AP (activity pattern; D-diurnal, C-cathemeral, N-nocturnal), CD (corneal diameter, mm), AL  (eye axial length, mm), Rel.Cornea (relative cornea size, CD/AL), Haplorhine (haplorhine primate, Y-yes, N-no)

- McDowell2024_UVA_data_subset.csv -- Subset of species from McDowell et al.'s lens transmission dataset, for calculating PropUVA and λT0.5.

- McDowell_UVA_T50_analysis.csv -- Estimated vs. published lens transmission metrics from McDowell et al. data, used in correlation analyses of species overlapping with Douglas and Jeffery's (2014) dataset. This is the output of Section 1 in the R code, with an Exclude column added ("New" = new taxa not in Douglas and Jeffery's data; Y = two species excluded from correlation (see Figure S1)).

- Hamster_digitized_data.csv -- Lens transmission curve digitized from Hut et al. (2000) Figure 5 for Syrian hamster (Mesocricetus auratus). Curve digitized using WebPlotDigitizer (version 5.2).

- mammal_trees_1000_Nov2026.nex -- subset of 1000 phylogenetic trees pruned to the 47 focal species that were used in phylogenetic regression analyses.

- Mammalian_Lens_Transmission.R -- Full analysis script: data preparation, phylogenetic regressions (Tables 1, 2, S4-S6), model diagnostics, and figure generation.


**Requirements:**
All analyses were run in R 4.6.1. Required packages: ape (v.5.8.1), car (v.3.1.5), dplyr (v.1.2.1), geiger (v.2.0.12), ggplot2 (v.4.0.3), phylolm (v.2.6.5), phytools (v.2.5.2), nationalparkcolors(v.0.1.0), but any color palette can be substituted.

**Tree provenance:**
The full set of 10,000 mammalian phylogenetic trees is from VertLife.org:
MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_FBDasZhouEtAl_all10k_v2_nexus.trees

From: Upham NS, Esselstyn JA, Jetz W. 2019. Inferring the mammal tree: Species-level sets of
phylogenies for questions in ecology, evolution, and conservation. PLOS Biology 17, e3000494.
(https://doi.org/10.1371/journal.pbio.3000494)

mammal_trees_1000_Nov2026.nex (included here) is a random 1000 tree subsample drawn from that set using set.seed(247) in R, then pruned to the 47 species in this dataset.

**Running the analysis:**
The R script is organized into six sections:
1. Calculate %UVA and 50Tnm -- calculates lens transmission metrics from McDowell et al. (2024) transmission curves and from the digitized hamster data.

2. Compare calculated transmission values with published values from Douglas & Jeffery (2014) for shared species/data to validate calculation methods.

3. Prepare Dataset and Trees -- loads and transforms the mammalian data in csv format. Loads in the 1000 trees.

4. Perform phylogenetic regressions -- runs the models reported in Tables 1-2, S4-S6. Expect console warnings during some analyses for models that drop taxa with missing data for the given formula. This is expected behavior because sample sizes differ by predictor and subset.

6. Model Diagnostic Checks -- checks variables in each model for collinearity using VIF, checks residuals of each model for normality and homoscedasticity, runs alternate models excluding specific taxa to confirm excluding them does not impact the results.

7. Making Figures -- generates Fig 1 panels, Fig. S1-S3.

**Contact:**
Carrie Veilleux (cveill AT midwestern.edu) for questions about the data or code.
