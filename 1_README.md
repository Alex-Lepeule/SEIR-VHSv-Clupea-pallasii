# SEIR-VHSv-Clupea-pallasii
A modified SEIR model to simulate the dynamic of the temperature dependent disease that is the viral hemorrhagic septicemia virus in pacific herring.

To start and estimate the parameters using a maximum likelihood model open the files names 'param estimates'. to match the parameters with the experimental data you will need to charge the .csv file named 'cumulative mortality' in this file you can find the value for the average cumulative mortality for each of the temperature tested. You can also find the standard deviation for each temperature. The results are from Hershberger et al. 2013. In the file named 'fitted parameters maximum likelihood' is the curves which represents the values of each parameter for the variation of temperature. 
To estimate the carrier induced infection rate, we had to do sensitivity analysis to aproximate the value of beta C, those data where not available and this was the best way to get close to the real value.

Once each parameters are done you can go the the modeling part. first you will need to open the file named 'SEIC_Constant_temperature.R' in the file is contained the simulation for each of the chanlenged temperature. 8.8°C, 11.2°C, and 14.7°C. Then you will have to open the 'Herd_immunity' file which contains the simulation with the proportion of carrier individulas as a function of the likelihood of outbreaks.

The file named 'SEIC_Temp_Variation.R' contains the model to simulate the impact of temperature variation. 

Then you can vary the number of carrier individulas to show the herd immunity in the case of temperature variation by oppening the 'SEIC_Tempe_and_Carrier-Variation.R'

The article in which those simulation have been used will be published soon, you can however use these codes for your work freely. 

Enjoy
