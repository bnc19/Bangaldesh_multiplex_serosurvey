Code and data to reproduce the analyses in:

Cracknell Daniels B, Kishor K.P, Vanhomwegen J, et al. The public health value of simultaneous assessment of serologic responses to multiple pathogens: Bangladesh as a case study.

The pipeline uses de-identified serosurvey data from Bangladesh (`data/deidentified_data.csv`), in which ages are grouped and coordinates rounded to two decimal places. Results produced from this code and data will therefore differ from those reported in the main manuscript.

The scripts are numbered and but can be run in any order to product estimates of the spatial clustering (script 1), risk-factors (scripts 2), and required study sample sizes (script 3).

Requires INLA: install.packages("INLA", repos=c(getOption("repos"), INLA="<https://inla.r-inla-download.org/R/stable>"), dep=TRUE)
