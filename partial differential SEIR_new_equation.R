
################################################################################
#############################                  #################################
######################## Temperature variation SEIC-P ##########################
#############################                  #################################
################################################################################

################################################################################
################################################################################
#############################  Required Library  ###############################
################################################################################
################################################################################

library(deSolve)
library(ggplot2)
library(truncnorm)
library(gridExtra)
library(grid)
library(lattice)
library(stats)
library(scales)
library(ggpubr)
library(patchwork)
library(plyr)
library(dplyr)
library(tidyr)

################################################################################
############################## Initial values ##################################
################################################################################

N = 100
E = 0
I = 1
C = 0
P = 0
D = 0
S=N-C-I-E-D

################################################################################
############################ Time step setting #################################
################################################################################

endT = 365
deltaT = 1
time = seq(0,endT,by = deltaT)
temp = 2.95*sin(time*(0.55/30))+11.75
################################################################################
######################## Parameters variation ##################################
################################################################################

beta <-matrix(0.0184,nrow = length(time))
sigma <- matrix(0.04199*temp-0.34936,nrow = length(time))
alpha <- matrix(-0.024049*temp+0.429996,nrow = length(time))
gamma <- matrix(0.058687*temp-0.453518,nrow = length(time))
tau <- matrix(0.01895*temp-0.04254,nrow = length(time))
delta <- matrix(-0.3959*temp+11.0789,nrow = length(time))
recoinf <- matrix(0.0046, nrow = length(time))
pie <- data.frame(recoinf, temp)
pie$recoinf <- ifelse(pie$temp > 9, 0, pie$recoinf)
pie <- pie$recoinf
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
                    pie,
                    phi,
                    mu,
                    rho)
################################################################################
########################### Event matrix setting ###############################
################################################################################

emat2 = matrix(c(-1,0,0,0,0,1,################################################## mu*S
                 -1,1,0,0,-1,0,################################################# beta * P * S
                 -1,1,0,0,-1,0,################################################# beta * P * S
                 0,-1,0,0,0,1,################################################## mu*E
                 0,-1,1,0,0,0,################################################## sigma*E
                 0,0,-1,0,0,1,################################################## mu*I
                 0,0,-1,0,0,1,################################################## alpha * I
                 0,0,-1,1,0,0,################################################## gamma*I
                 0,0,0,-1,0,1,################################################## mu*C
                 0,0,0,0,1,0,################################################### delta*I
                 0,0,0,0,1,0,################################################### phi*C
                 0,0,0,0,-1,0,################################################## tau*P
                 0,0,0,0,-1,0,################################################# rho*P
                 0,0,0,0,-1,0),##################################################-beta * P * S*10
               ncol = 6, byrow = TRUE)########################################## every possibility the model can take

################################################################################
########################## Equation setting ####################################
################################################################################

rlist2 = c(quote(ifelse(mu*S > S, S, mu*S)),
           quote(ifelse(beta*S*I > S, S, beta*S*I)),
           quote(ifelse(pie*S*C > S, S, pie*S*C)),
           quote(ifelse(mu*E > E, E, mu*E)),
           quote(ifelse(sigma*E > E, E, sigma*E)),
           quote(ifelse(mu*I > I, I, mu*I)),
           quote(ifelse(alpha * I > I, I, alpha * I)),
           quote(ifelse(gamma*I > I, I, gamma*I)),
           quote(ifelse(mu*C > C, C, mu*C)),
           quote(delta*I),
           quote(phi*C),
           quote(tau*P),
           quote(rho*P),
           quote(beta*S*P))

################################################################################
############################# Model function ###################################
################################################################################

SEICP = function(rateqs,eventmatrix, parameters, initialvals, deltaT, endT){
  time = seq(0,endT,by = deltaT)
  res <- data.frame(matrix(NA, ncol = length(initialvals)+ 1, nrow = length(time)))
  res[, 1] = time
  names(res) = c("time", names(inits))
  res[1, ] = c(0, inits)
  for(i in 1:(length(time)-1)){
    parameters <- c(beta = param$beta[i], sigma = param$sigma[i], gamma = param$gamma[i],
                    alpha = param$alpha[i], mu = param$mu[i], delta = param$delta[i],pie = param$pie[i],
                    phi = param$phi[i], tau = param$tau[i], rho = param$rho[i])
    #calculate overall rates
    rat = sapply(rateqs, eval, as.list(c(parameters, res[i,])))
    evts = rpois(1,sum(rat)*deltaT)
    if(evts > 0){
      #draw event
      whichevent = sample(1:nrow(eventmatrix), evts, prob = rat, replace = TRUE)
      mt = rbind(eventmatrix[whichevent,], t(matrix(res[i,-1])))
      mt = matrix(as.numeric(mt), ncol = ncol(mt))
      #update states
      res[i + 1, -1] = apply(mt, 2, sum)
      res[i + 1, ][res[i + 1,] < 0] = 0 
    }
    else{
      #if no events in deltaT
      res[i+1,-1] = res[i,-1]
    }
  }
  return(res)
}
################################################################################ Source: (Bjørnstad 2023)

################################################################################
#############################   var temperature  ###############################
################################################################################

inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
simTT <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)####################### apply the function
simTT <- data.frame(simTT)#################################################### creation of the data frame

simTT$iteration =0############################################################## Creation of the iterations
iterations = 499################################################################ set to 199 iterations more

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  simTT1 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  simTT1<-data.frame(simTT1)
  simTT1$iteration<-j
  simTT<-rbind(simTT1, simTT)
}############################################################################### for loop for iteration
simTT$iteration = as.factor(simTT$iteration)#################################### Set values for iteration as factor

################################################################################
################################################################################
#################################### Plots #####################################
################################################################################
################################################################################

SusTT <- data.frame(simTT$time,simTT$S)
transformed_SusTT <- SusTT %>%
  pivot_longer(cols = c("simTT.S"))
meansusTT <- ddply(transformed_SusTT, .variable = c("time","name"), summarise, values= mean(value))
simTT <- data.frame(simTT, meansusTT$values)

carTT <- data.frame(simTT$time,simTT$C)
transformed_carTT <- carTT %>%
  pivot_longer(cols = c("simTT.C"))
meancarTT <- ddply(transformed_carTT, .variable = c("time","name"), summarise, values= mean(value))
simTT <- data.frame(simTT, meancarTT$values)

expTT <- data.frame(simTT$time,simTT$E)
transformed_expTT <- expTT %>%
  pivot_longer(cols = c("simTT.E"))
meanexpTT <- ddply(transformed_expTT, .variable = c("time","name"), summarise, values= mean(value))
simTT <- data.frame(simTT, meanexpTT$values)

infTT <- data.frame(simTT$time,simTT$I)
transformed_infTT <- infTT %>%
  pivot_longer(cols = c("simTT.I"))
meaninfTT <- ddply(transformed_infTT, .variable = c("time","name"), summarise, values= mean(value))
simTT <- data.frame(simTT, meaninfTT$values)

parTT <- data.frame(simTT$time,simTT$P)
transformed_parTT <- parTT %>%
  pivot_longer(cols = c("simTT.P"))
meanparTT <- ddply(transformed_parTT, .variable = c("time","name"), summarise, values= mean(value))
simTT <- data.frame(simTT, meanparTT$values)
################################################################################
###################### Plot of the whole model #################################
################################################################################

TTplot <- ggplot(simTT,aes(x=time, by=iteration))+
  geom_line(aes(y=S),size = 3, colour = "blue4", alpha = 0.2)+
  geom_line(aes(y=E),size = 3, colour ="darkgoldenrod1", alpha = 0.2)+
  geom_line(aes(y=I),size = 3, colour = "brown1", alpha = 0.2)+
  geom_line(aes(y=C),size = 3, colour = "cyan3", alpha = 0.2)+
  geom_line(aes(y=D),size = 3, colour = "grey", alpha = 0.2)+
  ylab(label="Number of fish")+
  scale_y_continuous(labels = percent_format(scale = 1),
                     limits = c(0,N))+
  theme_classic()+
  xlab(label="Days post exposure")+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"),)

################################################################################
####################### Susceptible and carrier plot ###########################
################################################################################

SusceptibleTT <- ggplot(simTT,aes(x=time, by = iteration))+
  geom_line(aes(y=S, colour = "Susceptible"),size = 3, alpha = 0.2 )+
  geom_line(aes(y=meansusTT.values), colour = "darkblue", linewidth = 5)+
  ylab(label="
  Susceptible fish
  
       ")+
  scale_y_continuous(
    limits = c(0,N))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  xlab(label="Days post exposure")+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"))+
  theme(legend.justification=c(1,0), legend.position=c(0.6,0.5))+
  theme(legend.title=element_text(size=0,face="bold"),
        legend.background = element_rect(fill='#FFFFFF',
                                         size=100,linetype="solid"),
        legend.text=element_text(size=80),
        legend.key=element_rect(colour="#FFFFFF",
                                fill='#FFFFFF',
                                size=0,
                                linetype="solid"))+
  scale_colour_manual("",
                      breaks=c("Susceptible"),
                      values=c("blue4"))

CarrierTT <- ggplot(simTT,aes(x=time, by = iteration))+
  geom_line(aes(y=C, colour = "Carrier"),size = 3, alpha = 0.2 )+
  geom_line(aes(y=meancarTT.values), colour = "darkcyan", linewidth = 5)+
  ylab(label="
       carrier fish
       
       ")+
  scale_y_continuous(
    limits = c(0,N))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  xlab(label="Days post exposure")+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"))+
  theme(legend.justification=c(1,0), legend.position=c(0.6,0.5))+
  theme(legend.title=element_text(size=0,face="bold"),
        legend.background = element_rect(fill='#FFFFFF',
                                         size=100,linetype="solid"),
        legend.text=element_text(size=80),
        legend.key=element_rect(colour="#FFFFFF",
                                fill='#FFFFFF',
                                size=0,
                                linetype="solid"))+
  scale_colour_manual("",
                      breaks=c("Carrier"),
                      values=c("cyan3" ))
################################################################################
########################## Exposed and infected plot ###########################
################################################################################

InfectedTT <- ggplot(simTT,aes(x=time, by=iteration))+
  geom_line(aes(y=E, colour = "Exposed"),size = 3, alpha = 0.2)+
  geom_line(aes(y=I, colour = "Infected"),size = 3, alpha = 0.2)+
  geom_line(aes(y=meanexpTT.values), colour = "darkorange", linewidth = 5)+
  geom_line(aes(y=meaninfTT.values), colour = "darkred", linewidth = 5)+
  ylab(label="
  Exposed and 
       infected fish
       ")+
  scale_y_continuous(
    limits = c(0,N))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  xlab(label="Days post exposure")+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"))+
  theme(legend.justification=c(1,0), legend.position=c(0.5,0.6))+
  theme(legend.title=element_text(size=0,face="bold"),
        legend.background = element_rect(fill='#FFFFFF',
                                         size=100,linetype="solid"),
        legend.text=element_text(size=80),
        legend.key=element_rect(colour="#FFFFFF",
                                fill='#FFFFFF',
                                size=0,
                                linetype="solid"))+
  scale_colour_manual("",
                      breaks=c("Exposed","Infected"),
                      values=c("darkgoldenrod1", "brown1"))
################################################################################
########################## Exposed and infected plot ###########################
################################################################################

CarrierInfTT <- ggplot(simTT,aes(x=time, by = iteration))+
  geom_line(aes(y=C, colour = "Carrier"),size = 3, alpha = 0.2 )+
  geom_line(aes(y=I, colour = "Infected"),size = 3, alpha = 0.2)+
  geom_line(aes(y=meancarTT.values), colour = "darkcyan", linewidth = 5)+
  geom_line(aes(y=meaninfTT.values), colour = "darkred", linewidth = 5)+
  ylab(label="
  Infected and
       carrier fish
       ")+
  scale_y_continuous(
    limits = c(0,N))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  xlab(label="Days post exposure")+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"))+
  theme(legend.justification=c(1,0), legend.position=c(0.6,0.5))+
  theme(legend.title=element_text(size=0,face="bold"),
        legend.background = element_rect(fill='#FFFFFF',
                                         size=100,linetype="solid"),
        legend.text=element_text(size=80),
        legend.key=element_rect(colour="#FFFFFF",
                                fill='#FFFFFF',
                                size=0,
                                linetype="solid"))+
  scale_colour_manual("",
                      breaks=c("Carrier","Infected"),
                      values=c("cyan3","brown1" ))

################################################################################
############################## Particle plot ###################################
################################################################################

ParticlesTT <- ggplot(simTT,aes(x=time, by=iteration))+
  geom_line(aes(y=P),size = 3, colour = "grey27", alpha = 0.2)+
  geom_line(aes(y=meanparTT.values), colour = "black", size = 5)+
  ylab(label="
  Millions of Particles 
       in the water (in PFU/L)
       ")+
  xlab(label="Days post exposure")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"))

################################################################################
######################## Temperature variation plot ############################
################################################################################

temp <- data.frame(temp, time)
TT <- ggplot(temp, aes(x = time))+
  geom_line(aes(y = temp), colour = "black", size = 10)+
  ylab(label = "
  Temperature (°C)
  
  
       ")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  xlab(label="Days post exposure")+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"),
        title =  element_text(size = 80))

################################################################################
########################### creation of the plots ##############################
################################################################################

png(filename = "TimeTempPlot.png", height = 4500, width = 6400, res = 100 )##### name and resolution of the file

(TT| SusceptibleTT)/
  (CarrierTT|InfectedTT)############################################################ making a grid with the multiple plots
dev.off()####################################################################### send the plot on the chosen working directory

png(filename = "ParticlesXInfeced.png", height = 4500, width = 3200, res = 100 )##### name and resolution of the file

(CarrierInfTT)/
  (ParticlesTT)############################################################ making a grid with the multiple plots
dev.off()####################################################################### send the plot on the chosen working directory

