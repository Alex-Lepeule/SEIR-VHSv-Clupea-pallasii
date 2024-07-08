
##############################
############    ##############
#########  SEIR-P  ###########
############    ##############
##############################
library(deSolve)
library(ggplot2)
library(truncnorm)
library(gridExtra)
library(grid)
library(lattice)
library(stats)
library(ggpubr)
library(scales)
library(plyr)
library(dplyr)
library(tidyr)
library(patchwork)
siz = 3                



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
endT = 365
deltaT = 1
time = seq(0,endT,by = deltaT)
temp = 2.95*sin(time*(0.55/30))+11.75

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
#############
### 20% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.20
P = 0
D = 0
S=N-C-I-E

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho,param$pie)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
TTsim20 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
TTsim20 <- data.frame(TTsim20)

TTsim20$iteration =0
iterations = 499

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  TTsim201 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  TTsim201<-data.frame(TTsim201)
  TTsim201$iteration<-j
  TTsim20<-rbind(TTsim201, TTsim20)
}

TTsim20$iteration = as.factor(TTsim20$iteration)

Sus20 <- data.frame(TTsim20$time,TTsim20$S)
transformed_Sus20 <- Sus20 %>%
  pivot_longer(cols = c("TTsim20.S"))
meansus20 <- ddply(transformed_Sus20, .variable = c("time","name"), summarise, values= mean(value))
TTsim20 <- data.frame(TTsim20, meansus20$values)

car20 <- data.frame(TTsim20$time,TTsim20$C)
transformed_car20 <- car20 %>%
  pivot_longer(cols = c("TTsim20.C"))
meancar20 <- ddply(transformed_car20, .variable = c("time","name"), summarise, values= mean(value))
TTsim20 <- data.frame(TTsim20, meancar20$values)

inf20 <- data.frame(TTsim20$time,TTsim20$I)
transformed_inf20 <- inf20 %>%
  pivot_longer(cols = c("TTsim20.I"))
meaninf20 <- ddply(transformed_inf20, .variable = c("time","name"), summarise, values= mean(value))
TTsim20 <- data.frame(TTsim20, meaninf20$values)

exp20 <- data.frame(TTsim20$time,TTsim20$E)
transformed_exp20 <- exp20 %>%
  pivot_longer(cols = c("TTsim20.E"))
meanexp20 <- ddply(transformed_exp20, .variable = c("time","name"), summarise, values= mean(value))
TTsim20 <- data.frame(TTsim20, meanexp20$values)

part20 <- data.frame(TTsim20$time,TTsim20$P)
transformed_part20 <- part20 %>%
  pivot_longer(cols = c("TTsim20.P"))
meanpart20 <- ddply(transformed_part20, .variable = c("time","name"), summarise, values= mean(value))
TTsim20 <- data.frame(TTsim20, meanpart20$values)

dead20 <- data.frame(TTsim20$time,TTsim20$D)
transformed_dead20 <- dead20 %>%
  pivot_longer(cols = c("TTsim20.D"))
meandead20 <- ddply(transformed_dead20, .variable = c("time","name"), summarise, values= mean(value))
TTsim20 <- data.frame(TTsim20, meandead20$values)

################################################################
################  50%  #########################################
################################################################

N = 100
E = 0
I = 1
C = N*0.50
P = 0
D = 0
S=N-C-I-E

inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
TTsim50 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
TTsim50 <- data.frame(TTsim50)

TTsim50$iteration =0
iterations = 499

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  TTsim501 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  TTsim501<-data.frame(TTsim501)
  TTsim501$iteration<-j
  TTsim50<-rbind(TTsim501, TTsim50)
}
TTsim50$iteration = as.factor(TTsim50$iteration)

Sus50 <- data.frame(TTsim50$time,TTsim50$S)
transformed_Sus50 <- Sus50 %>%
  pivot_longer(cols = c("TTsim50.S"))
meansus50 <- ddply(transformed_Sus50, .variable = c("time","name"), summarise, values= mean(value))
TTsim50 <- data.frame(TTsim50, meansus50$values)

car50 <- data.frame(TTsim50$time,TTsim50$C)
transformed_car50 <- car50 %>%
  pivot_longer(cols = c("TTsim50.C"))
meancar50 <- ddply(transformed_car50, .variable = c("time","name"), summarise, values= mean(value))
TTsim50 <- data.frame(TTsim50, meancar50$values)

inf50 <- data.frame(TTsim50$time,TTsim50$I)
transformed_inf50 <- inf50 %>%
  pivot_longer(cols = c("TTsim50.I"))
meaninf50 <- ddply(transformed_inf50, .variable = c("time","name"), summarise, values= mean(value))
TTsim50 <- data.frame(TTsim50, meaninf50$values)

exp50 <- data.frame(TTsim50$time,TTsim50$E)
transformed_exp50 <- exp50 %>%
  pivot_longer(cols = c("TTsim50.E"))
meanexp50 <- ddply(transformed_exp50, .variable = c("time","name"), summarise, values= mean(value))
TTsim50 <- data.frame(TTsim50, meanexp50$values)

part50 <- data.frame(TTsim50$time,TTsim50$P)
transformed_part50 <- part50 %>%
  pivot_longer(cols = c("TTsim50.P"))
meanpart50 <- ddply(transformed_part50, .variable = c("time","name"), summarise, values= mean(value))
TTsim50 <- data.frame(TTsim50, meanpart50$values)

dead50 <- data.frame(TTsim50$time,TTsim50$D)
transformed_dead50 <- dead50 %>%
  pivot_longer(cols = c("TTsim50.D"))
meandead50 <- ddply(transformed_dead50, .variable = c("time","name"), summarise, values= mean(value))
TTsim50 <- data.frame(TTsim50, meandead50$values)

N = 100
E = 0
I = 1
C = N*0.80
P = 0
D = 0
S=N-C-I-E

inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
TTsim80 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
TTsim80 <- data.frame(TTsim80)

TTsim80$iteration =0
iterations = 499

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  TTsim801 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  TTsim801<-data.frame(TTsim801)
  TTsim801$iteration<-j
  TTsim80<-rbind(TTsim801, TTsim80)
}
TTsim80$iteration = as.factor(TTsim80$iteration)

Sus80 <- data.frame(TTsim80$time,TTsim80$S)
transformed_Sus80 <- Sus80 %>%
  pivot_longer(cols = c("TTsim80.S"))
meansus80 <- ddply(transformed_Sus80, .variable = c("time","name"), summarise, values= mean(value))
TTsim80 <- data.frame(TTsim80, meansus80$values)

car80 <- data.frame(TTsim80$time,TTsim80$C)
transformed_car80 <- car80 %>%
  pivot_longer(cols = c("TTsim80.C"))
meancar80 <- ddply(transformed_car80, .variable = c("time","name"), summarise, values= mean(value))
TTsim80 <- data.frame(TTsim80, meancar80$values)

inf80 <- data.frame(TTsim80$time,TTsim80$I)
transformed_inf80 <- inf80 %>%
  pivot_longer(cols = c("TTsim80.I"))
meaninf80 <- ddply(transformed_inf80, .variable = c("time","name"), summarise, values= mean(value))
TTsim80 <- data.frame(TTsim80, meaninf80$values)

exp80 <- data.frame(TTsim80$time,TTsim80$E)
transformed_exp80 <- exp80 %>%
  pivot_longer(cols = c("TTsim80.E"))
meanexp80 <- ddply(transformed_exp80, .variable = c("time","name"), summarise, values= mean(value))
TTsim80 <- data.frame(TTsim80, meanexp80$values)

part80 <- data.frame(TTsim80$time,TTsim80$P)
transformed_part80 <- part80 %>%
  pivot_longer(cols = c("TTsim80.P"))
meanpart80 <- ddply(transformed_part80, .variable = c("time","name"), summarise, values= mean(value))
TTsim80 <- data.frame(TTsim80, meanpart80$values)

dead80 <- data.frame(TTsim80$time,TTsim80$D)
transformed_dead80 <- dead80 %>%
  pivot_longer(cols = c("TTsim80.D"))
meandead80 <- ddply(transformed_dead80, .variable = c("time","name"), summarise, values= mean(value))
TTsim80 <- data.frame(TTsim80, meandead80$values)




N = 100
E = 0
I = 1
C = N*0.90
P = 0
D = 0
S=N-C-I-E

inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
TTsim90 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
TTsim90 <- data.frame(TTsim90)

TTsim90$iteration =0
iterations = 499

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  TTsim901 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  TTsim901<-data.frame(TTsim901)
  TTsim901$iteration<-j
  TTsim90<-rbind(TTsim901, TTsim90)
}
TTsim90$iteration = as.factor(TTsim90$iteration)

Sus90 <- data.frame(TTsim90$time,TTsim90$S)
transformed_Sus90 <- Sus90 %>%
  pivot_longer(cols = c("TTsim90.S"))
meansus90 <- ddply(transformed_Sus90, .variable = c("time","name"), summarise, values= mean(value))
TTsim90 <- data.frame(TTsim90, meansus90$values)

car90 <- data.frame(TTsim90$time,TTsim90$C)
transformed_car90 <- car90 %>%
  pivot_longer(cols = c("TTsim90.C"))
meancar90 <- ddply(transformed_car90, .variable = c("time","name"), summarise, values= mean(value))
TTsim90 <- data.frame(TTsim90, meancar90$values)

inf90 <- data.frame(TTsim90$time,TTsim90$I)
transformed_inf90 <- inf90 %>%
  pivot_longer(cols = c("TTsim90.I"))
meaninf90 <- ddply(transformed_inf90, .variable = c("time","name"), summarise, values= mean(value))
TTsim90 <- data.frame(TTsim90, meaninf90$values)

exp90 <- data.frame(TTsim90$time,TTsim90$E)
transformed_exp90 <- exp90 %>%
  pivot_longer(cols = c("TTsim90.E"))
meanexp90 <- ddply(transformed_exp90, .variable = c("time","name"), summarise, values= mean(value))
TTsim90 <- data.frame(TTsim90, meanexp90$values)

part90 <- data.frame(TTsim90$time,TTsim90$P)
transformed_part90 <- part90 %>%
  pivot_longer(cols = c("TTsim90.P"))
meanpart90 <- ddply(transformed_part90, .variable = c("time","name"), summarise, values= mean(value))
TTsim90 <- data.frame(TTsim90, meanpart90$values)

dead90 <- data.frame(TTsim90$time,TTsim90$D)
transformed_dead90 <- dead90 %>%
  pivot_longer(cols = c("TTsim90.D"))
meandead90 <- ddply(transformed_dead90, .variable = c("time","name"), summarise, values= mean(value))
TTsim90 <- data.frame(TTsim90, meandead90$values)

#######################################
##Plots 20
#######################################

SusceptibleTT20 <- ggplot(TTsim20,aes(x=time, by = iteration))+
  geom_line(aes(y=S),size = 3, colour = "blue3", alpha = 0.2)+
  geom_line(aes(y=C),size = 3, colour = "cyan", alpha = 0.2)+
  geom_line(aes(y=meansus20.values), size = 5, colour = "darkblue")+
  geom_line(aes(y=meancar20.values), size = 5, colour = "darkcyan")+
  ylab(label="Susceptible and 
  carrier fish
       ")+
  ggtitle("Temperature variation 
          with 20% of carriers")+
  scale_y_continuous(
    limits = c(0,N))+
  theme_classic()+   
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  xlab(label="")+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"),
        title =  element_text(size = 95, face = "bold"),
        axis.title.y = element_text(angle = 0))

InfectedTT20 <- ggplot(TTsim20,aes(x=time, by=iteration))+
  geom_line(aes(y=E),size = siz, colour = "darkgoldenrod1", alpha = 0.2)+
  geom_line(aes(y=I),size = siz, colour = "brown1", alpha = 0.2)+
  geom_line(aes(y=meanexp20.values), size = 5, colour = "darkorange2")+
  geom_line(aes(y=meaninf20.values), size = 5, colour = "darkred")+
  ylab(label=" Infected and 
  exposed fish
       ")+
  scale_y_continuous(
    limits = c(0,N))+
  theme_classic()+   theme(panel.border = element_rect(color = "black", fill = NA, size = 6),         panel.background = element_rect(fill = "white"))+
  xlab(label="")+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"),
        axis.title.y = element_text(angle = 0))

DeadTT20 <- ggplot(TTsim20,aes(x=time, by=iteration))+
  geom_line(aes(y=D),size = siz, colour = "grey", alpha = 0.2)+
  geom_line(aes(y=meandead20.values), size = 5, colour = "black")+
  ylab(label="Cumulative
       mortality")+
  scale_y_continuous(
    limits = c(0,N))+
  theme_classic()+   theme(panel.border = element_rect(color = "black", fill = NA, size = 6),         panel.background = element_rect(fill = "white"))+
  xlab(label="")+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"),
        axis.title.y = element_text(angle = 0))

ParticlesTT20 <- ggplot(TTsim20,aes(x=time, by=iteration))+
  geom_line(aes(y=P),size = 3, colour = "grey27", alpha = 0.2)+
  geom_line(aes(y=meanpart20.values), size = 5, colour = "black")+
  ylab(label="Total particles 
  shed in the water 
       (in 10^6 PFU)")+
  xlab(label="Days post exposure")+
  scale_y_continuous(limits = c(0,600))+
  theme_classic()+   theme(panel.border = element_rect(color = "black", fill = NA, size = 6),         panel.background = element_rect(fill = "white"))+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"),
        axis.title.y = element_text(angle = 0))



##################################################
################# Plots 50
################################################

SusceptibleTT50 <- ggplot(TTsim50,aes(x=time, by = iteration))+
  geom_line(aes(y=S),size = 3, colour = "blue4", alpha = 0.2)+
  geom_line(aes(y=C),size = 3, colour = "cyan3", alpha = 0.2)+
  geom_line(aes(y=meansus50.values), size = 5, colour = "darkblue")+
  geom_line(aes(y=meancar50.values), size = 5, colour = "darkcyan")+
  ylab(label="
       ")+
  ggtitle("Temperature variation 
          with 50% of carriers")+
  scale_y_continuous(
    limits = c(0,N))+
  theme_classic()+   theme(panel.border = element_rect(color = "black", fill = NA, size = 6),         panel.background = element_rect(fill = "white"))+
  xlab(label="")+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"),
        title =  element_text(size = 95, face = "bold"))

InfectedTT50 <- ggplot(TTsim50,aes(x=time, by=iteration))+
  geom_line(aes(y=E),size = siz, colour = "darkgoldenrod1", alpha = 0.2)+
  geom_line(aes(y=I),size = siz, colour = "brown1", alpha = 0.2)+
  geom_line(aes(y=meanexp50.values), size = 5, colour = "darkorange2")+
  geom_line(aes(y=meaninf50.values), size = 5, colour = "darkred")+
  ylab(label="
       ")+
  scale_y_continuous(
    limits = c(0,N))+
  theme_classic()+   theme(panel.border = element_rect(color = "black", fill = NA, size = 6),         panel.background = element_rect(fill = "white"))+
  xlab(label="")+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"))

DeadTT50 <- ggplot(TTsim50,aes(x=time, by=iteration))+
  geom_line(aes(y=D),size = siz, colour = "grey", alpha = 0.2)+
  geom_line(aes(y=meandead50.values), size = 5, colour = "black")+
  ylab(label="
       ")+
  scale_y_continuous(
    limits = c(0,N))+
  theme_classic()+   theme(panel.border = element_rect(color = "black", fill = NA, size = 6),         panel.background = element_rect(fill = "white"))+
  xlab(label="")+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"))

ParticlesTT50 <- ggplot(TTsim50,aes(x=time, by=iteration))+
  geom_line(aes(y=P),size = 3, colour = "grey27", alpha = 0.2)+
  geom_line(aes(y=meanpart50.values), size = 5, colour = "black")+
  ylab(label="
       ")+
  xlab(label="Days post exposure")+
  scale_y_continuous(limits = c(0,600))+
  theme_classic()+   theme(panel.border = element_rect(color = "black", fill = NA, size = 6),         panel.background = element_rect(fill = "white"))+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"))


###########################################################
##################### plots 80
###########################################################

SusceptibleTT80 <- ggplot(TTsim80,aes(x=time, by = iteration))+
  geom_line(aes(y=S),size = 3, colour = "blue4", alpha = 0.2)+
  geom_line(aes(y=C),size = 3, colour = "cyan3", alpha = 0.2)+
  geom_line(aes(y=meansus80.values), size = 5, colour = "darkblue")+
  geom_line(aes(y=meancar80.values), size = 5, colour = "darkcyan")+
  ylab(label="
       ")+
  ggtitle("Temperature variation 
          with 80% of carriers")+
  scale_y_continuous(
    limits = c(0,N))+
  theme_classic()+   theme(panel.border = element_rect(color = "black", fill = NA, size = 6),         panel.background = element_rect(fill = "white"))+
  xlab(label="")+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=75,face="bold"),
        title =  element_text(size = 95, face = "bold"))

InfectedTT80 <- ggplot(TTsim80,aes(x=time, by=iteration))+
  geom_line(aes(y=E),size = siz, colour = "darkgoldenrod1", alpha = 0.2)+
  geom_line(aes(y=I),size = siz, colour = "brown1", alpha = 0.2)+
  geom_line(aes(y=meanexp80.values), size = 5, colour = "darkorange2")+
  geom_line(aes(y=meaninf80.values), size = 5, colour = "darkred")+
  ylab(label="
       ")+
  scale_y_continuous(
    limits = c(0,N))+
  theme_classic()+   theme(panel.border = element_rect(color = "black", fill = NA, size = 6),         panel.background = element_rect(fill = "white"))+
  xlab(label="")+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"))

DeadTT80 <- ggplot(TTsim80,aes(x=time, by=iteration))+
  geom_line(aes(y=D),size = siz, colour = "grey", alpha = 0.2)+
  geom_line(aes(y=meandead80.values), size = 5, colour = "black")+
  ylab(label="
       ")+
  scale_y_continuous(
    limits = c(0,N))+
  theme_classic()+   theme(panel.border = element_rect(color = "black", fill = NA, size = 6),         panel.background = element_rect(fill = "white"))+
  xlab(label="")+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"))

ParticlesTT80 <- ggplot(TTsim80,aes(x=time, by=iteration))+
  geom_line(aes(y=P),size = 3, colour = "grey27", alpha = 0.2)+
  geom_line(aes(y=meanpart80.values), size = 5, colour = "black")+
  ylab(label="
       ")+
  xlab(label="Days post exposure")+
  scale_y_continuous(limits = c(0,600))+
  theme_classic()+   theme(panel.border = element_rect(color = "black", fill = NA, size = 6),         panel.background = element_rect(fill = "white"))+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"))

#################################################
######### Plot 90
################################################

SusceptibleTT90 <- ggplot(TTsim90,aes(x=time, by = iteration))+
  geom_line(aes(y=S),size = 3, colour = "blue4", alpha = 0.2)+
  geom_line(aes(y=C),size = 3, colour = "cyan3", alpha = 0.2)+
  geom_line(aes(y=meansus90.values), size = 5, colour = "darkblue")+
  geom_line(aes(y=meancar90.values), size = 5, colour = "darkcyan")+
  ylab(label="
       ")+
  ggtitle("Temperature variation 
          with 90% of carriers")+
  scale_y_continuous(
    limits = c(0,N))+
  theme_classic()+   theme(panel.border = element_rect(color = "black", fill = NA, size = 6),         panel.background = element_rect(fill = "white"))+
  xlab(label="")+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=75,face="bold"),
        title =  element_text(size = 95, face = "bold"))

InfectedTT90 <- ggplot(TTsim90,aes(x=time, by=iteration))+
  geom_line(aes(y=E),size = siz, colour = "darkgoldenrod1", alpha = 0.2)+
  geom_line(aes(y=I),size = siz, colour = "brown1", alpha = 0.2)+
  geom_line(aes(y=meanexp90.values), size = 5, colour = "darkorange2")+
  geom_line(aes(y=meaninf90.values), size = 5, colour = "darkred")+
  ylab(label="
       ")+
  scale_y_continuous(
    limits = c(0,N))+
  theme_classic()+   theme(panel.border = element_rect(color = "black", fill = NA, size = 6),         panel.background = element_rect(fill = "white"))+
  xlab(label="")+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"))


DeadTT90 <- ggplot(TTsim90,aes(x=time, by=iteration))+
  geom_line(aes(y=D),size = siz, colour = "grey", alpha = 0.2)+
  geom_line(aes(y=meandead90.values), size = 5, colour = "black")+
  ylab(label="
       ")+
  scale_y_continuous(
    limits = c(0,N))+
  theme_classic()+   theme(panel.border = element_rect(color = "black", fill = NA, size = 6),         panel.background = element_rect(fill = "white"))+
  xlab(label="")+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"))

ParticlesTT90 <- ggplot(TTsim90,aes(x=time, by=iteration))+
  geom_line(aes(y=P),size = 3, colour = "grey27", alpha = 0.2)+
  geom_line(aes(y=meanpart90.values), size = 5, colour = "black")+
  ylab(label="
       ")+
  xlab(label="Days post exposure")+
  
  scale_y_continuous(limits = c(0,600))+
  theme_classic()+   theme(panel.border = element_rect(color = "black", fill = NA, size = 6),         panel.background = element_rect(fill = "white"))+
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100,face="bold"))


png(filename = "carrier 20 50 80 90.png", height = 9000, width = 12800, res = 100 )

(SusceptibleTT20|SusceptibleTT50|SusceptibleTT80|SusceptibleTT90)/
  (InfectedTT20|InfectedTT50|InfectedTT80|InfectedTT90)/
  (DeadTT20|DeadTT50|DeadTT80|DeadTT90)

dev.off()
