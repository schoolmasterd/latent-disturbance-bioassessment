
This code accompanies the manuscript:
A new method for bioassassement for ecosystems with complex communities and environmental gradients

Donald R. Schoolmaster Jr.<sup>1</sup>, Valerie Partridge<sup>2</sup>

1 U.S. Geological Survey, Wetland and Aquatic Research Center. email: schoolmasterd@usgs.gov

2 State of Washington, Department of Ecology

The Code folder contains four files.
The script "calculate_d_scores.R" follows the manuscript's section "Case study: development of a Puget Sound Marine Benthic Index." It uses the designated training data to calculate D scores for the training period and estimate and save the parameters necessary for updating with other data. This script contains code for creating figures 2-4 of the manuscript

The script "adjust_d_for_stressors.R" follows the manuscript's subsection "Using D to identify and quantify stressors." It shows a worked example of adjusting the D scores for TOC and TN to determine what proportion of the D scores can be attributed to those stressors. This script contains code for creating Figures 6 and 7 and Table 2.

The script "metric_comparisons.R" calculates metics and creates figure 5

The script "Fig_makers.R" contains functions for plotting results used by the other scripts.  


The Data folder contains 3 files and a folder.
The file "habitat and TOC and species abundances for Long-term 2017 and 2018 samples.csv" is the code used for estimating the parameters for the model.

The file "sediment chemistry for 2017-2018 baseline.csv" contains the data on stressors.

The file "AMBI EG for 551 taxa in 2017-2018.csv" contains EG categorizations for the taxa for calculating an AMBI for figure 5

The folder "Model" is the destination of the parameters needed for updating the models. These files are created by the script "calculate_d_scores.R"

The Output folder is the main destination for the figures created by the scripts.

USGS IP-155119

approx. runtime: 10 min
