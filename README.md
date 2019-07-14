# Retro-Pink

## Data RP4
Experimental data is available as a single csv file (sep ','). The data is pre-processed: trials with RTs < 150ms were removed, as well as all trials displaying unexpected timings and other signs of hardware/software failure.\
The file contains a header indicating what each column corresponds to.\
1st column indicates participant ID\
2nd column contains response times (ms). NaN values indicate trials where participants opted-out.\
3rd column is trial accuracy. 1 if correct (hit or correct rejection), 0 if incorrect (miss or false alarm)\
4th column is trial congruence. 1 if cue and target were on the same side, 0 if not. In cue-absent or target-absent trials, congruence cannot be assessed and a NaN value is recorded. \
5th column is SOA in ms. As for congruence, NaN values were recorded when less than two stimuli were presented.\
6th column indicates target location on-screen. 1 if Left, 2 if Right, NaN if target-absent.\
7th column indicates cue location. 1 if Left, 2 if Right, NaN if cue-absent.\

## Matlab .m script
This script was used to produce Figure 1 and classical statistical analyses performed outside modeling. The Statistics toolbox may be necessary to run this script. In theory, all that is needed to run it is to provide the appropriate path to the RP4 dataset.
Note that a small bug may occur when trying to create the figures. In this case, execute the section line by line.

## R .R script
This script is an adaptation of the script from Anders et al. 2016 as indicated in the body of the article. It was used to produce Figure 2 and Supplementary Figure 1. In theory, you should only need to provide the appropriate data path to the RP4 dataset to use this script. On a recent machine, the whole process takes under 30 minutes to run.

## Contact
Please contact Daphné Rimsky-Robert if you have questions, suggestions or feedback at drr@tuta.io.
