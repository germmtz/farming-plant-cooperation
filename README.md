# farming-plant-cooperation
Code used for the study "Farming plant cooperation in crops"

Two approaches were used in the study: analytical predictions and individual-based simulations

The analytical approach is developed in the notebook "Farming_plant_cooperation_analytical_model.nb". We used Wolfram Mathematica v. 11.3.0.0

Individual-based simulations were performed with R software v. 3.5.2. The main code is available in the file "Farming_plant_cooperation_ind_based_simulations.R". This file uses "Farming_plant_cooperation_fecundity.R" and "Farming_plant_cooperation_R_converted_analytical_model.R" as source code

"Farming_plant_cooperation_fecundity.R" contains the function that compute plant fecundity

"Farming_plant_cooperation_R_converted_analytical_model.R" contains functions derived from the analytical approach (relatedness, selection gradient, etc). These functions are used to compare the simulation results with the theoretical predictions
