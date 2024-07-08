
####Library required ################

library(ggplot2)
library(gridExtra)
library(grid)
library(scales)
library(patchwork)

#############################################################################################################################
######################## Set of the time step and the variation of temperature as function of time ##########################
#############################################################################################################################

endT = 365 #length
deltaT = 1 #time step
time = seq(0,endT,by = deltaT)
temp = 2.95*sin(time*(0.55/30))+11.75 # temp = f(time)

#############################################################################################################################
############################################## Parameters variation #########################################################
#############################################################################################################################

beta <-matrix(0.0184,nrow = length(time))
sigma <- matrix(0.04199*temp-0.34936,nrow = length(time))
alpha <- matrix(-0.024049*temp+0.429996,nrow = length(time))
gamma <- matrix(0.058687*temp-0.453518,nrow = length(time))
tau <- matrix(0.01895*temp-0.04254,nrow = length(time))
delta <- matrix(-0.3959*temp+11.0789,nrow = length(time))
recoinf <- matrix(0.0046, nrow = length(time))
pie <- data.frame(recoinf, temp)
pie$recoinf <- ifelse(pie$temp > 9, 0, pie$recoinf)
pie <- pie$recoshed
recoshed <- matrix(3.5, nrow = length(time))
f <- data.frame(recoshed, temp)
f$recoshed <- ifelse(f$temp > 9, 0, f$recoshed)
phi <- f$recoshed
mu<- matrix(0.0008, ncol= 1, nrow = length(time))
rho<- matrix(0.9, ncol= 1, nrow = length(time))
param <- data.frame(beta, 
                    sigma, 
                    alpha, 
                    gamma,
                    tau,
                    delta,
                    phi,
                    mu,
                    rho)
summary(lm(c(0.5,1,1.5)~c(8.8,11.2,14.7)))

#############################################################################################################################
############################################### Mortality rate ##############################################################
#############################################################################################################################

# Creation of the data frame
sigmaE <- data.frame(sigma, temp) 
###### ggplot setting
sigmaPlot <- ggplot(sigmaE, aes(x = temp))+
  # set the value of alpha for each constant temperature
  # the parameters are calculated from Hershberger et al. 2013 
  geom_line(aes(y = sigma), colour = "brown1", size = 10)+
  geom_point(aes(y = 0.6543,x = 8.8), colour = "black", size =30)+ 
  geom_point(aes(y = 0.9284,x = 11.2), colour = "black", size =30)+ 
  geom_point(aes(y = 0.9223,x = 14.7), colour = "black", size =30)+
  ######## Set name of Y axis
  ylab(label = "
  σ (fish*day-1)
       ")+
  ############ panel style
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 6),
        panel.background = element_rect(fill = "white"))+
  theme(axis.line = element_line(linewidth = 2))+
  ######## Set name of X axis
  xlab(label = "
  Temperature (°C)
       ")+
  #scale on high definition
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"),
        title =  element_text(size = 100))

# linear regression and fitted Rsquare
summary(lm(c(0.6543,0.9284,0.9223)~c(8.8,11.2,14.7)))

#############################################################################################################################
############################################### Mortality rate ##############################################################
#############################################################################################################################

# Creation of the data frame
alphaE <- data.frame(alpha, temp) 
###### ggplot setting
alphaPlot <- ggplot(alphaE, aes(x = temp))+
  # set the value of alpha for each constant temperature
  # the parameters are calculated from Hershberger et al. 2013 
  geom_line(aes(y = alpha), colour = "brown1", size = 10)+
  geom_point(aes(y = 0.2060,x = 8.8), colour = "black", size =30)+ 
  geom_point(aes(y = 0.1815,x = 11.2), colour = "black", size =30)+ 
  geom_point(aes(y = 0.0680,x = 14.7), colour = "black", size =30)+
  ######## Set name of Y axis
  ylab(label = "
  α (fish*day-1)
       ")+
  ############ panel style
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 6),
        panel.background = element_rect(fill = "white"))+
  theme(axis.line = element_line(linewidth = 2))+
  ######## Set name of X axis
  xlab(label = "
  Temperature (°C)
       ")+
  #scale on high definition
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"),
        title =  element_text(size = 100))

# linear regression and fitted Rsquare
summary(lm(c(0.2060,0.1815,0.0680)~c(8.8,11.2,14.7)))

#############################################################################################################################
############################################### Recovery rate ###############################################################
#############################################################################################################################

# Creation of the data frame
gammaE <- data.frame(gamma, temp) 
###### ggplot setting
gammaPlot <- ggplot(gammaE, aes(x = temp))+
  # set the value of gamma for each constant temperature
  # the parameters are calculated from Hershberger et al. 2013 
  geom_line(aes(y = gamma), colour = "brown1", size = 10)+
  geom_point(aes(y = 0.0472,x = 8.8), colour = "black", size =30)+ 
  geom_point(aes(y = 0.2303,x = 11.2), colour = "black", size =30)+ 
  geom_point(aes(y = 0.3984,x = 14.7), colour = "black", size =30)+
  ######## Set name of Y axis
  ylab(label = "
  γ (fish*day-1)
       ")+
  ############ panel style
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth =  6),
        panel.background = element_rect(fill = "white"))+                   
  theme(axis.line = element_line(linewidth = 2))+
  ######## Set name of X axis
  xlab(label = "
  Temperature (°C)
       ")+
  #scale on high definition
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"),
        title =  element_text(size = 100))
# linear regression and fitted Rsquare
summary(lm(c(0.0472,0.2303,0.3984)~c(8.8,11.2,14.7)))

#############################################################################################################################
############################################### Infected shedding rate ######################################################
#############################################################################################################################

# Creation of the data frame
cE <- data.frame(delta, temp)
###### ggplot setting
cPlot <- ggplot(cE, aes(x = temp))+
  # set the value of teta for each constant temperature
  # the parameters are calculated from Hershberger et al. 2013
  geom_line(aes(y = delta), colour = "brown1", size = 10)+
  geom_point(aes(y = 7.8,x = 8.8), colour = "black", size =30)+
  geom_point(aes(y = 6.3,x = 11.2), colour = "black", size =30)+ 
  geom_point(aes(y = 5.4,x = 14.7), colour = "black", size =30)+
  ######## Set name of Y axis
  ylab(label = "
  δ (10^6PFU/fish/day)
       ")+
  ############ panel style
  theme_classic() +
  scale_y_continuous(labels = label_number(accuracy = 0.01),
                     limits = c(5,8))+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 6),
        panel.background = element_rect(fill = "white"))+
  theme(axis.line = element_line(linewidth = 2))+
  ######## Set name of X axis
  xlab(label = "
       Temperature (°C)
       ")+
  #scale on high definition
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"),
        title =  element_text(size = 100))
# linear regression and fitted Rsquare
summary(lm(c(7.8,6.3,5.4)~c(8.8,11.2,14.7)))

#############################################################################################################################
############################################### Inactivation rate ###########################################################
#############################################################################################################################

# Creation of the data frame
rE <- data.frame(tau, temp)
###### ggplot setting
rPlot <- ggplot(rE, aes(x = temp))+
  # set the value of tau for each constant temperature
  # the parameters are calculated from Hawley and Garver 2008
  geom_line(aes(y = tau), colour = "brown1", size = 10)+
  geom_point(aes(y = 0.13,x = 8.8), colour = "black", size =30)+
  geom_point(aes(y = 0.16,x = 11.2), colour = "black", size =30)+
  geom_point(aes(y = 0.24,x = 14.7), colour = "black", size =30)+
  ######## Set name of Y axis
  ylab(label = "
  τ (10^6PFU*day-1)
       ")+
  ############ panel style
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 6),
        panel.background = element_rect(fill = "white"))+
  theme(axis.line = element_line(linewidth = 2))+
  ######## Set name of X axis
  xlab(label = "
  Temperature (°C)
       ")+
  #scale on high definition
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"), 
        title =  element_text(size = 100))
# linear regression and fitted Rsquare
summary(lm(c(0.13,0.16,0.24)~c(8.8,11.2,14.7)))

#############################################################################################################################
############################################### Carrier shedding rate #######################################################
#############################################################################################################################

# Creation of the data frame
fE <- data.frame(phi,temp)
###### ggplot setting
fPlot <- ggplot(fE, aes(x = temp))+
  # set the value of phi for each constant temperature
  # the parameters are calculated from Hershberger et al. 2021
  geom_line(aes(y = phi), colour = "brown1", size = 10)+
  geom_point(aes(y = 1.5,x = 8.8), colour = "black", size =30)+
  geom_point(aes(y = 0,x = 11.2), colour = "black", size =30)+  
  geom_point(aes(y = 0,x = 14.7), colour = "black", size =30)+
  ######## Set name of Y axis
  ylab(label = "
  φ (10^6PFU/fish/day)
       ")+
  ############ panel style
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 6),
        panel.background = element_rect(fill = "white"))+
  theme(axis.line = element_line(linewidth = 2))+
  ######## Set name of X axis
  xlab(label = "
  Temperature (°C)
       ")+
  #scale on high definition
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"),
        title =  element_text(size = 100))


#############################################################################################################################
############################################### creation of the plots #######################################################
#############################################################################################################################

png(filename = "Param.png", height = 4500, width = 6400, res = 100 )# name and resolution of the file
# making a grid with the multiple plots
(sigmaPlot|alphaPlot)/
             (gammaPlot|fPlot)/
             (cPlot|rPlot) 
dev.off() #send the plot on the chosen working directory

#beta
0.0551
#sigma
summary(lm(c(0.6543,0.9284,0.9223)~c(8.8,11.2,14.7)))
#gamma
summary(lm(c(0.0472,0.2303,0.3984)~c(8.8,11.2,14.7)))
#alpha
summary(lm(c(0.2060,0.1815,0.0680)~c(8.8,11.2,14.7)))
