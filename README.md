# Distinct neural subspaces during obstacle avoidance

This repository contains the code for the analysis and generating the figures for the publication:

Xing D, Truccolo W, Borton DA. “Emergence of distinct neural subspaces in motor cortical dynamics during volitional adjustments of ongoing locomotion”. *J Neurosci*. 2022 Oct 25:JN-RM-0746-22. doi: 10.1523/JNEUROSCI.0746-22.2022.

The main code consists of the 9 *makeFig_…m* files, each corresponding to a figure in the paper:

- *makeFig_Setup.m* - Figure 2 panel c
- *makeFig_Kinematics.m* - Figure 2 panels d-i
- *makefig_NeuralDataSupp.m* - Figure 3
- *makeFig_exampleNeurons.m* - Figure 4
- *makeFig_NeuralProps.m* - Figure 5
- *makeFig_dPCA.m* - Figure 7
- *makeFig_PLDS.m* - Figure 8
- *makeFig_ModelPLDS.m* - Figure 9
- *makeFig_TimeCourse.m* - Figure 10

The data analyses carried out in the paper (e.g. dPCA, calculating dispersion, decoding, ect) are all implemented at the beginning of the corresponding *makeFig_…m* file. The exception is the analysis for fitting the Poisson Linear Dynamical Systems (PLDS) models to the data, which is carried out in *runPLDS.m* instead, and the analysis simulating the neural data from the kinematics, which is carried out in *addKinModelSpikes.m*

To run the code, the data containing the neural recordings and kinematics needs to be located in a directory ./Data. The datasets TrialsDataBoomer.mat and TrialsDataStarbuck.mat can be requested from [David_Borton@brown.edu](mailto:David_Borton@brown.edu)

The code utilizes previously published toolboxes for dPCA, jPCA, and circular statistics (referenced in the paper), which should be downloaded and added to the matlab path in order to carry out the analyses. 

The code for fitting the PLDS model is included in the PLDS_Model folder. Additionally, helper functions for some of the functions are required and can be obtained from [https://github.com/davidxingy/NeuralProcessingTools](https://github.com/davidxingy/NeuralProcessingTools)