###########################################################################
30 Mar 2022
David Lusseau

We use this fork to present the code alteration to DEPONS proposed in Lusseau, Kindt-Larsen and van Beest (2022) Emergent interactions in the management of multiple threats to the conservation objective of harbour porpoises [url to come]. 
as well as deposit the R code for the statistical modelling carried out on the simulation outcomes (new folder on this fork named "outcome_analyses"). 

The simulation outcomes for the scenarios described in the paper and used in the subsequent statistical analyses can be found at doi:10.11583/DTU.19455554

Alterations to the DEPONS code to introduce condition mediation of avoidance can be found in the folder "CM_porpoise". THe folder contains 3 versions of the porpoise.java file with altered avoidance response as described in equations 2-4 of the paper.


_Therefore_, this fork is not synced with the DEPONS master link.
###########################################################################


# DEPONS 2.2
Simulating effects of disturbances on harbour porpoises in the North Sea.

This is the second public version of the DEPONS model for simulation of population effects of noise for the harbour porpoise. It differs from version 1.1 (the first public version) in being calibrated entirely based on data from the North Sea (see TRACE document for details).

The objective of the DEPONS model is to simulate how harbor porpoise population dynamics are affected by pile-driving noise associated with construction of offshore wind farms. The animals’ survival is directly related to their energy levels, and the population dynamics are affected by noise through its impact on the animals’ foraging behavior. By ensuring that the animals’ movement patterns, space use and reactions to noise are realistic, the population dynamics in the model have the same causal drivers in the model as in nature.

The model was developed as part of the DEPONS project (Disturbance Effects of POrpoises in the North Sea) at Aarhus University, Denmark. The code was developed in Eclipse/Repast, and must be installed under this framework. See http://www.depons.dk or contact Jacob Nabe-Nielsen (jnn@bios.au.dk) for more information.
