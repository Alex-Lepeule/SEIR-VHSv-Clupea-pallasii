################################################################################
##########################                  ####################################
###############  SEIC-P for each constant temperature  #########################
##########################                  ####################################
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
library(ggpubr)
library(scales)
library(plyr)
library(dplyr)
library(tidyr)
library(patchwork)

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

endT = 25
deltaT = 1
time = seq(0,endT,by = deltaT)

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

rlist2 = c(quote(ifelse(mu*S > S, S, mu*S)),#natural mortality susceptible
           quote(ifelse(beta*S*I > S, S, beta*S*I)), #transition of susceptible to exposed initiated by infected
           quote(ifelse(phi*S*C > S, S, phi*S*C)), #transition of susceptible to exposed initiated by Carriers
           quote(ifelse(mu*E > E, E, mu*E)),#natural mortality Exposed
           quote(ifelse(sigma*E > E, E, sigma*E)),#transition from Exposed to infected (latency)
           quote(ifelse(mu*I > I, I, mu*I)),#Natural mortality infected 
           quote(ifelse(alpha * I > I, I, alpha * I)), #disease induced mortality
           quote(ifelse(gamma*I > I, I, gamma*I)),#Transition from Infected to carrier
           quote(ifelse(mu*C > C, C, mu*C)),#Natural mortality carrier
           quote(delta*I), #viral shedding of infected 
           quote(phi*C), #viral shedding of carrier
           quote(tau*P),#viral inactivation
           quote(rho*P),#viral flushing 
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
#############################   cold temperature  ##############################
################################################################################

inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
paras = c(mu = 0.0008,beta = 0.0184, sigma = 0.6543, gamma = 0.0472, 
         alpha = 0.2060, delta = 7.8, phi = 0.0046, tau = 0.13, rho = 0.9)######## Parameter values
simCold = SEICP(rlist2, emat2, paras, inits,deltaT,endT)######################## apply the function
simCold <- data.frame(simCold)################################################## creation of the data frame
#alpha/gamma
0.104/0.017

simCold$iteration = 0 ########################################################## Creation of the iterations
iterations = 499 ############################################################### set to 199 iterations more

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  simCold1 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  simCold1<-data.frame(simCold1)
  simCold1$iteration<-j
  simCold<-rbind(simCold1, simCold)
} ############################################################################## for loop for iteration
simCold$iteration = as.factor(simCold$iteration)################################ Set values for iteration as factor





################################################################################
########################### Ambient temperature ################################
################################################################################

inits = c(S =S, E =E, I =I, C =C , P = P, D = D)################################ Initial values

paras = c(mu = 0.0008, beta =  0.0184, sigma = 0.9284, gamma = 0.2303, 
       alpha = 0.1815, delta = 6.3, phi = 0,tau = 0.16, rho = 0.9)########### Parameter values
simAmbient = SEICP(rlist2, emat2, paras, inits,deltaT,endT)##################### apply the function
simAmbient <- data.frame(simAmbient)############################################ creation of the data frame


simAmbient$iteration = 0 ####################################################### Creation of the iterations
iterations = 499 ############################################################### set to 199 iterations more

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  simAmb1 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  simAmb1<-data.frame(simAmb1)
  simAmb1$iteration<-j
  simAmbient<-rbind(simAmb1, simAmbient)
} ############################################################################## for loop for iteration
simAmbient$iteration = as.factor(simAmbient$iteration)########################## Set values for iteration as factor


################################################################################
############################  warm temperature #################################
################################################################################

inits = c(S =S, E =E, I =I, C =C , P = P, D = D) ############################### Initial values

paras = c(mu = 0.0008, beta =  0.0184, sigma = 0.9223, gamma = 0.3984, 
        alpha = 0.068, delta = 5.4, phi = 0, tau = 0.24, rho = 0.9)######### Parameter values
simWarm = SEICP(rlist2, emat2, paras, inits,deltaT,endT)######################## apply the function
simWarm <- data.frame(simWarm)################################################## creation of the data frame


simWarm$iteration = 0 ########################################################## Creation of the iterations
iterations = 499 ############################################################### set to 199 iterations more

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  simWarm1 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  simWarm1<-data.frame(simWarm1)
  simWarm1$iteration<-j
  test <- subset(simWarm1, simWarm1$time == endT)
  ifelse(simWarm1$C < 100,simWarm1,simWarm1 == matrix(NA, ncol = ncol(simWarm1), nrow = nrow(simWarm1)))
  simWarm<-rbind(simWarm1, simWarm)
} ############################################################################## for loop for iteration
simWarm$iteration = as.factor(simWarm$iteration)################################ Set values for iteration as factor

################################################################################
################################################################################
#################################### Plots #####################################
################################################################################
################################################################################

all <- data.frame(simWarm,simAmbient, simCold)################################## set data frame with all models results

################################################################################
############################ Susceptible plot ##################################
################################################################################
Suscold <- data.frame(all$time,all$S.2)
transformed_cold <- Suscold %>%
  pivot_longer(cols = c("all.S.2"))
meancold <- ddply(transformed_cold, .variable = c("time","name"), summarise, values= mean(value))

all <- data.frame(all, meancold$values)
allplotSCold <- ggplot(all,aes(x=time,by = c(iteration)))+
  geom_line(aes(y=S.2),size = 2, colour = "blue", alpha = 0.2)+
  geom_line(aes(y=meancold.values),colour = "darkblue",size = 5)+
  ylab(label="
  Susceptible
  fish
  ")+
  scale_y_continuous(limits = c(0,N))+
  xlab(label="")+
  ggtitle("Cold temperature
        ~ 8.8°C")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  theme(axis.text=element_text(size=50),
        axis.title=element_text(size=75,face="bold"),
        axis.title.y = element_text(angle = 0),
        title = element_text(size = 70, face = "bold"))+
  theme(legend.justification=c(1,0), legend.position=c(0.9,0.5))+
  theme(legend.title=element_text(size=30,face="bold"),
        legend.background = element_rect(fill='#FFFFFF',
                                         size=0.5,linetype="solid"),
        legend.text=element_text(size=80),
        legend.key=element_rect(colour="#FFFFFF",
                                fill='#FFFFFF',
                                linewidth = 1,
                                linetype="solid"))

SusAmbient <- data.frame(all$time,all$S.1)
transformed_amb <- SusAmbient %>%
  pivot_longer(cols = c("all.S.1"))
meanAmbient <- ddply(transformed_amb, .variable = c("time","name"), summarise, values= mean(value))
all <- data.frame(all, meanAmbient$values)

allplotSAmbient <- ggplot(all,aes(x=time,by = c(iteration)))+
  geom_line(aes(y=S.1),colour = "orange",size = 2, alpha = 0.2)+
  geom_line(aes(y=meanAmbient.values),colour = "darkorange4",size = 5)+
  ylab(label="")+
  scale_y_continuous(limits = c(0,N))+
  xlab(label="")+
  ggtitle("Ambient temperature
        ~ 11.2°C")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  theme(axis.text=element_text(size=50),
        axis.title=element_text(size=50,face="bold"),
        title = element_text(size = 70, face = "bold"))+
  theme(legend.justification=c(1,0), legend.position=c(0.9,0.5))+
  theme(legend.title=element_text(size=30,face="bold"),
        legend.background = element_rect(fill='#FFFFFF',
                                         size=0.5,linetype="solid"),
        legend.text=element_text(size=80),
        legend.key=element_rect(colour="#FFFFFF",
                                fill='#FFFFFF',
                                linewidth = 1,
                                linetype="solid"))

SusWarm <- data.frame(all$time,all$S)
transformed_warm <- SusWarm %>%
  pivot_longer(cols = c("all.S"))
meanWarm <- ddply(transformed_warm, .variable = c("time","name"), summarise, values= mean(value))
all <- data.frame(all, meanWarm$values)

allplotSWarm <- ggplot(all,aes(x=time,by = c(iteration)))+
  geom_line(aes(y=S),colour = "red",size = 2, alpha = 0.2)+
  geom_line(aes(y=meanWarm.values),colour = "darkred",size = 5)+
  ylab(label="
       ")+
  scale_y_continuous(limits = c(0,N))+
  xlab(label="")+
  ggtitle("Warm temperature
        ~ 14.7°C")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  theme(axis.text=element_text(size=50),
        axis.title=element_text(size=50,face="bold"),
        title = element_text(size = 70, face = "bold"))+
  theme(legend.justification=c(1,0), legend.position=c(0.9,0.5))+
  theme(legend.title=element_text(size=30,face="bold"),
        legend.background = element_rect(fill='#FFFFFF',
                                         size=0.5,linetype="solid"),
        legend.text=element_text(size=80),
        legend.key=element_rect(colour="#FFFFFF",
                                fill='#FFFFFF',
                                linewidth = 1,
                                linetype="solid"))

################################################################################
############################## Infected plot ###################################
################################################################################

Infcold <- data.frame(all$time,all$I.2)
transformedI_cold <- Infcold %>%
  pivot_longer(cols = c("all.I.2"))
meanIcold <- ddply(transformedI_cold, .variable = c("time","name"), summarise, values= mean(value))
all <- data.frame(all, meanIcold$values)

allplotICold <- ggplot(all,aes(x=time,by = c(iteration)))+
  geom_line(aes(y=I.2),size = 2, colour = "blue", alpha = 0.2)+
  geom_line(aes(y=meanIcold.values),colour = "darkblue",size = 5)+
  ylab(label=" Infected
  fish
       ")+
  xlab(label="")+
  scale_y_continuous(limits = c(0,N))+
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  theme(axis.text=element_text(size=50),
        axis.title=element_text(size=75,face="bold"),
        axis.title.y = element_text(angle = 0))

InfAmbient <- data.frame(all$time,all$I.1)
transformedI_ambient <- InfAmbient %>%
  pivot_longer(cols = c("all.I.1"))
meanIAmbient <- ddply(transformedI_ambient, .variable = c("time","name"), summarise, values= mean(value))
all <- data.frame(all, meanIAmbient$values)

allplotIAmbient <- ggplot(all,aes(x=time,by = c(iteration)))+
  geom_line(aes(y=I.1),size = 2, colour = "orange", alpha = 0.2)+
  geom_line(aes(y=meanIAmbient.values),colour = "darkorange4",size = 5)+
  ylab(label="")+
  xlab(label="")+
  scale_y_continuous(limits = c(0,N))+
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  theme(axis.text=element_text(size=50),
        axis.title=element_text(size=50,face="bold"))

InfWarm <- data.frame(all$time,all$I)
transformedI_warm <- InfWarm %>%
  pivot_longer(cols = c("all.I"))
meanIWarm <- ddply(transformedI_warm, .variable = c("time","name"), summarise, values= mean(value))
all <- data.frame(all, meanIWarm$values)

allplotIWarm <- ggplot(all,aes(x=time,by = c(iteration)))+
  geom_line(aes(y=I),size = 2, colour = "red", alpha = 0.2)+
  geom_line(aes(y=meanIWarm.values),colour = "darkred",size = 5)+
  ylab(label="")+
  xlab(label="")+
  scale_y_continuous(limits = c(0,N))+
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  theme(axis.text=element_text(size=50),
        axis.title=element_text(size=50,face="bold"))

################################################################################
############################### Carrier plot ###################################
################################################################################

Carcold <- data.frame(all$time,all$C.2)
transformedC_cold <- Carcold %>%
  pivot_longer(cols = c("all.C.2"))
meanCcold <- ddply(transformedC_cold, .variable = c("time","name"), summarise, values= mean(value))

all <- data.frame(all, meanCcold$values)

allplotCCold <- ggplot(all,aes(x=time,by = iteration))+
  geom_line(aes(y=C.2),size = 2, colour = "blue", alpha = 0.2)+
  geom_line(aes(y=meanCcold.values),colour = "darkblue",size = 5)+
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  ylab(label=" Carrier
  fish
       ")+
  xlab(label="")+
  scale_y_continuous(limits = c(0,N))+
  theme(axis.text=element_text(size=50),
        axis.title=element_text(size=75,face="bold"),
        axis.title.y = element_text(angle = 0))

CarAmbient <- data.frame(all$time,all$C.1)
transformedC_amb <- CarAmbient %>%
  pivot_longer(cols = c("all.C.1"))
meanCAmbient <- ddply(transformedC_amb, .variable = c("time","name"), summarise, values= mean(value))
all <- data.frame(all, meanCAmbient$values)

allplotCAmbient <- ggplot(all,aes(x=time,by = iteration))+
  geom_line(aes(y=C.1),size = 2, colour = "orange", alpha = 0.2)+
  geom_line(aes(y=meanCAmbient.values),colour = "darkorange4",size = 5)+
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  ylab(label="")+
  xlab(label="")+
  scale_y_continuous(limits = c(0,N))+
  theme(axis.text=element_text(size=50),
        axis.title=element_text(size=50,face="bold"))

CarWarm <- data.frame(all$time,all$C)
transformedC_warm <- CarWarm %>%
  pivot_longer(cols = c("all.C"))
meanCWarm <- ddply(transformedC_warm, .variable = c("time","name"), summarise, values= mean(value))
all <- data.frame(all, meanCWarm$values)

allplotCWarm <- ggplot(all,aes(x=time,by = iteration))+
  geom_line(aes(y=C),size = 2, colour = "red", alpha = 0.2)+
  geom_line(aes(y=meanCWarm.values),colour = "darkred",size = 5)+
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  ylab(label="")+
  xlab(label="")+
  scale_y_continuous(limits = c(0,N))+
  theme(axis.text=element_text(size=50),
        axis.title=element_text(size=50,face="bold"))

################################################################################
############################### Mortality plot #################################
################################################################################

Deathcold <- data.frame(all$time,all$D.2)
transformedD_cold <- Deathcold %>%
  pivot_longer(cols = c("all.D.2"))
meanDcold <- ddply(transformedD_cold, .variable = c("time","name"), summarise, values= mean(value))
all <- data.frame(all, meanDcold$values)

allplotDCold <- ggplot(all,aes(x=time,by = iteration))+
  geom_line(aes(y=D.2),size = 2, colour = "blue", alpha = 0.2)+
  geom_line(aes(y=meanDcold.values),colour = "darkblue",size = 5)+
  ylab(label="Cumulative 
  mortality
       ")+
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  xlab(label="")+
  scale_y_continuous(limits = c(0,N))+
  theme(axis.text=element_text(size=50),
        axis.title=element_text(size=75,face="bold"),
        axis.title.y = element_text(angle = 0))

DeathAmbient <- data.frame(all$time,all$D.1)
transformedD_amb <- DeathAmbient %>%
  pivot_longer(cols = c("all.D.1"))
meanDAmbient <- ddply(transformedD_amb, .variable = c("time","name"), summarise, values= mean(value))
all <- data.frame(all, meanDAmbient$values)

allplotDAmbient <- ggplot(all,aes(x=time,by = iteration))+
  geom_line(aes(y=D.1),size = 2, colour = "orange", alpha = 0.2)+
  geom_line(aes(y=meanDAmbient.values),colour = "darkorange4",size = 5)+
  ylab(label="")+
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  xlab(label="")+
  scale_y_continuous(limits = c(0,N))+
  theme(axis.text=element_text(size=50),
        axis.title=element_text(size=50,face="bold"))

DeathWarm <- data.frame(all$time,all$D)
transformedD_warm <- DeathWarm %>%
  pivot_longer(cols = c("all.D"))
meanDWarm <- ddply(transformedD_warm, .variable = c("time","name"), summarise, values= mean(value))
all <- data.frame(all, meanDWarm$values)

allplotDWarm <- ggplot(all,aes(x=time,by = iteration))+
  geom_line(aes(y=D),size = 2, colour = "red", alpha = 0.2)+
  geom_line(aes(y=meanDWarm.values),colour = "darkred",size = 5)+
  ylab(label="")+
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  xlab(label="")+
  scale_y_continuous(limits = c(0,N))+
  theme(axis.text=element_text(size=50),
        axis.title=element_text(size=50,face="bold"))

################################################################################
############################## Particle plot ###################################
################################################################################
Partcold <- data.frame(all$time,all$P.2)
transformedP_cold <- Partcold %>%
  pivot_longer(cols = c("all.P.2"))
meanPcold <- ddply(transformedP_cold, .variable = c("time","name"), summarise, values= mean(value))
all <- data.frame(all, meanPcold$values)

allplotPCold <- ggplot(all,aes(x=time,by = iteration))+
  geom_line(aes(y=P.2),size = 2, colour = "blue", alpha = 0.2)+
  geom_line(aes(y=meanPcold.values),colour = "darkblue",size = 5)+
  ylab(label="Total particles
  shed in the water 
  (10^6 PFUs)
       ")+
  theme_classic() +
  scale_y_continuous(limits = c(0,500))+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  xlab(label="Days post exposure")+
  theme(axis.text=element_text(size=50),
        axis.title=element_text(size=75,face="bold"),
        axis.title.y = element_text(angle = 0))

PartAmbient <- data.frame(all$time,all$P.1)
transformedP_amb <- PartAmbient %>%
  pivot_longer(cols = c("all.P.1"))
meanPAmbient <- ddply(transformedP_amb, .variable = c("time","name"), summarise, values= mean(value))
all <- data.frame(all, meanPAmbient$values)

allplotPAmbient <- ggplot(all,aes(x=time,by = iteration))+
  geom_line(aes(y=P.1),size = 2, colour = "orange", alpha = 0.2)+
  geom_line(aes(y=meanPAmbient.values),colour = "darkorange4",size = 5)+
  ylab(label="")+
  theme_classic() +
  scale_y_continuous(limits = c(0,500))+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  xlab(label="Days post exposure")+
  theme(axis.text=element_text(size=50),
        axis.title=element_text(size=75,face="bold"))

PartWarm <- data.frame(all$time,all$P)
transformedP_warm <- PartWarm %>%
  pivot_longer(cols = c("all.P"))
meanPWarm <- ddply(transformedP_warm, .variable = c("time","name"), summarise, values= mean(value))
all <- data.frame(all, meanPWarm$values)

allplotPWarm <- ggplot(all,aes(x=time,by = iteration))+
  geom_line(aes(y=P),size = 2, colour = "red", alpha = 0.2)+
  geom_line(aes(y=meanPWarm.values),colour = "darkred",size = 5)+
  ylab(label="")+
  theme_classic() +
  scale_y_continuous(limits = c(0,500))+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  xlab(label="Days post exposure")+
  theme(axis.text=element_text(size=50),
        axis.title=element_text(size=75,face="bold"))

textCold <- text_grob("
Cold simulation
8.8°C
          ",size = 100, face = "bold")
textAmbient <- text_grob("
Ambient simulation
11.2°C
          ",size = 100, face = "bold")
textwarm <- text_grob("
Warm simulation
14.7°C
          ",size = 100, face = "bold")

textSus <- text_grob("
Susceptible fish

          ",size = 100, face = "bold")
textinf <- text_grob("
Infected fish

       ",size = 100, face = "bold")
textCar <- text_grob("
Carrier fish

          ",size = 100, face = "bold")
textDead <- text_grob("
Dead fish

       ",size = 100, face = "bold")
textPar <- text_grob("Total particles 
  shed in the water 
       (in 10^6 PFU)",size = 100, face = "bold")
textwhite <- text_grob("")
################################################################################
########################### creation of the plots ##############################
################################################################################

png(filename = "AllPlot2.png", height = 4500, width = 6400, res = 100 )########## name and resolution of the file

(allplotSCold|allplotSAmbient|allplotSWarm)/
  (allplotICold|allplotIAmbient|allplotIWarm)/
  (allplotCCold|allplotCAmbient|allplotCWarm)/
  (allplotDCold|allplotDAmbient|allplotDWarm)
############################################################ making a grid with the multiple plots
dev.off()####################################################################### send the plot on the chosen working directory

