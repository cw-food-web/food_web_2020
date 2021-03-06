# food_web_2020
This repository holds scripts, datasets, and an Rmarkdown for the 2020 lake trout food web project. The contents of this repository are named and briefly described below.

The Rmarkdown provides a step-by-step explanation of the analyses in this project.\
rmd_food_web_10012020_14.html (Rmarkdown html document explaining unified R script's analyses with code chunks and figures)\
rmd_food_web_10012020_14.rmd (the R file to produce the html version of the Rmarkdown)

The unified R script is annotated, calls the study's main isotope dataset, and includes code to reproduce analyses.\
unified_food_web_script_08242020_29.R (study's main R script)\
isotope_dataset_MT_lake_trout_food_web_11072020.csv (study's main dataset)

Besides the main isotope dataset, the unified R script also calls three other datasets:\
transect_data.csv  (for NMDS ordinations)\
lake_data.csv  (for NMDS ordinations)\
conversion_09302020_2.csv (for binomial regression)

The lake trout diet MixSIAR script is separate from the unified script.\
lake_trout_mixsiar_model_08272020_3.R (lake trout diet R script)\
lake_trout_mix.csv (MixSIAR mixtures; this is called in the lake trout diet R script)\
lake_trout_source2.csv (MixSIAR sources; this is called in the lake trout diet R script)\
lake_trout_trophic_discrimination2.csv (MixSIAR discrimination factors; this is called in the lake trout diet R script)

Two other scripts make multi-panel plots from the contents of the unified food web script.\
trophic_disruption_superplot_11132020.R (combines trophic dispersion and displacement results into one figure)\
script_superplot_overlap_conversion_11132020_3.R (combines diet overlap and conversion results into one figure)
