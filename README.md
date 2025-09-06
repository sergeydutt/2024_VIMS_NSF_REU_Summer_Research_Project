Overview:
This repository contains code and materials from my NSF REU marine science research internship at the Virginia Institute of Marine Science. This project focused on calibrating a phytoplankton primary production model for the Chesapeake Bay and evaluating its performance across seasons and regions. Our hypothesis was that various seasons and regions in the Chesapeake Bay come with unique phytoplankton characteristics that require their own set of model assumptions. 

Repository Contents:
Paper_Figuers.R - script that was used to create all of the figures present in the paper and presentation
train_data_VA_AP_MD_all_06102024.RData - train data set for calibrating the model. This data set covers the Chesapeake Bay from 1984-2009 and has over 13,000 data points.
Sergey_Dutt_FinalPresentation.pdf - A pdf version of my final research presentation. This is a good summary of my research project.
Dutt_Sergey_FinalPaper.pdf - A pdf of my final research paper. This paper covers my summer research project more thoroughly and in depth than the presentation. 

Requirements:
At least R version: 4.4.0 (tested on macOS)
Packages: install.packages(c("rsample", "purrr", "infer", "Metrics", "tidyverse", "ggplot2”))

Usage Statement:
This repository showcases the work I completed during my NSF REU internship at VIMS. You are welcome to explore the paper and presentation, and to run the R script with the provided dataset to reproduce the figures. The script can be modified to test different model skill metrics or to experiment with calibrating additional models. I also welcome suggestions for more efficient coding practices—feel free to email me with ideas or revisions, as long as the script continues to achieve its original goals.

To run script: source("code/analysis.R")
