Data and code related to the article “Modeling temperature-dependant herd immunity and survival of Pacific herring exposed to viral haemorrhagic septicemia virus” and submitted in Ecosphere 2025

Data:

1_description_project.txt: Background explanation of the study and description of the model and how to exploit files and codes linked in the github repository.

2_cumulative mortality.csv: file containing the cumulative mortality data from Hershebrger et al 2013 (Hershberger, P. K., M. K. Purcell, L. M. Hart, et al. 2013. ‘Influence of Temperature on Viral Hemorrhagic Septicemia (Genogroup IVa) in Pacific Herring, Clupea Pallasii Valenciennes’. Journal of Experimental Marine Biology and Ecology 444 (June): 81–86. https://doi.org/10.1016/j.jembe.2013.03.006.) used as validation data for the maximum likelihood model. 
Column: left to right, for the three temperatures for which data were available. 
mean8C = the average value of cumulative mortality at each time step evidenced in Hershberger et al. 2013
Lower SD 8C = lower limit of the standard deviation for the challenges at 8°C
Upper SD 8C = Upper limit of the standard deviation for the challenges at 8°C.
Each group of three follows the same logic for the other challenged temperatures 11°C and 14°C  
This file is loaded in the 3_param_estimates.R

Code:

3_param_estimates.R: in this code we used the validation data from Hershberger et al. 2013 to best fit the model's dynamic for each temperature to the observed data using a maximum likelihood model.

4_fited_parameters_max_likelihood.R: Graphical representation of the variation of the  temperature dependent parameters after fitting to validation data.

5_SEIC_constant_temperature.R: First model with stochastic parameters and specific parameters for each simulated constant temperature. 

6_SEIC_Temp_Variation.R: In this script we simulated a temperature variation and varied the parameters according to the temperature. The simulations are starting in the spring at average temperature and go through a summer and winter period to simulate the full scale of our data.

7_SEIC_Temp_and_Carrier-Variation.R: In this script, we kept the temperature variation from the previous simulations and added a variation in the number of Recovered Convalesced Carriers (RCC) 4 different starting conditions are simulated 20%, 50%, 80% and 90% of RCCs at the start of the simulations.

8_Herd_immunity.R: In this script we varied the proportion of RCCs from 0% to 95% for each temperature to follow the evolution of the likelihood of herd immunity for each starting condition and each temperature.

9_Sensitivity_Analysis.R: Estimation of the value of BetaC which is the infection rate of RCC herring while being infectious under 9°C. We followed the variation of likelihood of herd immunity according to BetaC’ value. 
