
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

endT = 30
deltaT = 1
time = seq(0,endT,by = deltaT)
temp = matrix(9,nrow = length(time))

beta <-matrix(0.0184,nrow = length(time))
sigma <- matrix(0.04199*temp+0.34936,nrow = length(time))
alpha <- matrix(-0.024049*temp+0.429996,nrow = length(time))
gamma <- matrix(0.058687*temp-0.453518,nrow = length(time))
tau <- matrix(0.01895*temp-0.04254,nrow = length(time))
delta <- matrix(-0.3959*temp+11.0789,nrow = length(time))
recoinf <- matrix(0.0184, nrow = length(time))
pie <- data.frame(recoinf, temp)
pie$recoinf <- ifelse(pie$temp >= 9.1, 0, pie$recoinf)
pie <- pie$recoinf
recoshed <- matrix(3.5, nrow = length(time))
f <- data.frame(recoshed, temp)
f$recoshed <- ifelse(f$temp >= 9.1, 0, f$recoshed)
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
### 0% of recovered in the population
#############
N = 100
E = 0
I = 1
C = N*0.0
P = 0
D = 0
S=N-C-I-E

iter = 100
iterations = 99

coldsim02 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                        "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                        "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim0 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim0 <- data.frame(coldsim0)
  
  coldsim0$iteration =0
  
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim01 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim01<-data.frame(coldsim01)
    coldsim01$iteration<-j
    coldsim0<-rbind(coldsim01, coldsim0)
  }
  coldsim0$iter <- i
  coldsim02 <- rbind(coldsim02, coldsim0)
}


#############
### 10% of recovered in the population
#############
N = 100
E = 0
I = 1
C = N*0.10
P = 0
D = 0
S=N-C-I-E

iter = 100

coldsim102 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim10 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim10 <- data.frame(coldsim10)
  
  coldsim10$iteration =0
  
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim101 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim101<-data.frame(coldsim101)
    coldsim101$iteration<-j
    coldsim10<-rbind(coldsim101, coldsim10)
  }
  coldsim10$iter <- i
  coldsim102 <- rbind(coldsim102, coldsim10)
}


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
iter = 100

coldsim202 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim20 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim20 <- data.frame(coldsim20)
  
  coldsim20$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim201 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim201<-data.frame(coldsim201)
    coldsim201$iteration<-j
    coldsim20<-rbind(coldsim201, coldsim20)
  }
  coldsim20$iter <- i
  coldsim202 <- rbind(coldsim202, coldsim20)
}

#############
### 30% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.30
P = 0
D = 0
S=N-C-I-E

coldsim302 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim30 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim30 <- data.frame(coldsim30)
  
  coldsim30$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim301 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim301<-data.frame(coldsim301)
    coldsim301$iteration<-j
    coldsim30<-rbind(coldsim301, coldsim30)
  }
  coldsim30$iter <- i
  coldsim302 <- rbind(coldsim302, coldsim30)
}


#############
### 40% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.40
P = 0
D = 0
S=N-C-I-E

coldsim402 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim40 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim40 <- data.frame(coldsim40)
  
  coldsim40$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim401 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim401<-data.frame(coldsim401)
    coldsim401$iteration<-j
    coldsim40<-rbind(coldsim401, coldsim40)
  }
  coldsim40$iter <- i
  coldsim402 <- rbind(coldsim402, coldsim40)
}



#############
### 50% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.50
P = 0
D = 0
S=N-C-I-E

coldsim502 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim50 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim50 <- data.frame(coldsim50)
  
  coldsim50$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim501 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim501<-data.frame(coldsim501)
    coldsim501$iteration<-j
    coldsim50<-rbind(coldsim501, coldsim50)
  }
  coldsim50$iter <- i
  coldsim502 <- rbind(coldsim502, coldsim50)
}


#############
### 60% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.60
P = 0
D = 0
S=N-C-I-E

coldsim602 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim60 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim60 <- data.frame(coldsim60)
  
  coldsim60$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim601 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim601<-data.frame(coldsim601)
    coldsim601$iteration<-j
    coldsim60<-rbind(coldsim601, coldsim60)
  }
  coldsim60$iter <- i
  coldsim602 <- rbind(coldsim602, coldsim60)
}

#############
### 70% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.70
P = 0
D = 0
S=N-C-I-E

coldsim702 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim70 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim70 <- data.frame(coldsim70)
  
  coldsim70$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim701 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim701<-data.frame(coldsim701)
    coldsim701$iteration<-j
    coldsim70<-rbind(coldsim701, coldsim70)
  }
  coldsim70$iter <- i
  coldsim702 <- rbind(coldsim702, coldsim70)
}


#############
### 80% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.80
P = 0
D = 0
S=N-C-I-E

coldsim802 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim80 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim80 <- data.frame(coldsim80)
  
  coldsim80$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim801 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim801<-data.frame(coldsim801)
    coldsim801$iteration<-j
    coldsim80<-rbind(coldsim801, coldsim80)
  }
  coldsim80$iter <- i
  coldsim802 <- rbind(coldsim802, coldsim80)
}



#############
### 90% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.90
P = 0
D = 0
S=N-C-I-E

coldsim902 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim90 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim90 <- data.frame(coldsim90)
  
  coldsim90$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim901 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim901<-data.frame(coldsim901)
    coldsim901$iteration<-j
    coldsim90<-rbind(coldsim901, coldsim90)
  }
  coldsim90$iter <- i
  coldsim902 <- rbind(coldsim902, coldsim90)
}

cold0<- coldsim02 %>% filter(coldsim02$time == 30)
cold0 <- cold0 %>%
  filter(S > 1) %>%
  count(iteration)

cold10<- coldsim102 %>% filter(coldsim102$time == 30)
cold10 <- cold10 %>%
  filter(S > 1) %>%
  count(iteration)

cold20<- coldsim202 %>% filter(coldsim202$time == 30)
cold20 <- cold20 %>%
  filter(S > 1) %>%
  count(iteration)

cold30<- coldsim302 %>% filter(coldsim302$time == 30)
cold30 <- cold30 %>%
  filter(S > 1) %>%
  count(iteration)

cold40<- coldsim402 %>% filter(coldsim402$time == 30)
cold40 <- cold40 %>%
  filter(S > 1) %>%
  count(iteration)

cold50<- coldsim502 %>% filter(coldsim502$time == 30)
cold50 <- cold50 %>%
  filter(S > 1) %>%
  count(iteration)

cold60<- coldsim602 %>% filter(coldsim602$time == 30)
cold60 <- cold60 %>%
  filter(S > 1) %>%
  count(iteration)

cold70<- coldsim702 %>% filter(coldsim702$time == 30)
cold70 <- cold70 %>%
  filter(S > 1) %>%
  count(iteration)

cold80<- coldsim802 %>% filter(coldsim802$time == 30)
cold80 <- cold80 %>%
  filter(S > 1) %>%
  count(iteration)

cold90<- coldsim902 %>% filter(coldsim902$time == 3)
cold90 <- cold90 %>%
  filter(S > 1) %>%
  count(iteration)


cold0SD <-sd(cold0$n)
cold0MEAN <-mean(cold0$n)
coldCI0 <- 1.5/sqrt(100) 

cold10SD <-sd(cold10$n)
cold10MEAN <-mean(cold10$n)
coldCI10 <- cold10SD/sqrt(100)

cold20SD <-sd(cold20$n)
cold20MEAN <-mean(cold20$n)
coldCI20 <- cold20SD/sqrt(100) 

cold30SD <-sd(cold30$n)
cold30MEAN <-mean(cold30$n)
coldCI30 <- cold30SD/sqrt(100) 

cold40SD <-sd(cold40$n)
cold40MEAN <-mean(cold40$n)
coldCI40 <- cold40SD/sqrt(100) 

cold50SD <-sd(cold50$n)
cold50MEAN <-mean(cold50$n)
coldCI50 <- cold50SD/sqrt(100) 

cold60SD <-sd(cold60$n)
cold60MEAN <-mean(cold60$n)
coldCI60 <- cold60SD/sqrt(100) 

cold70SD <-sd(cold70$n)
cold70MEAN <-mean(cold70$n)
coldCI70 <- cold70SD/sqrt(100) 

cold80SD <-sd(cold80$n)
cold80MEAN <-mean(cold80$n)
coldCI80 <- cold80SD/sqrt(100) 

cold90SD <-sd(cold90$n)
cold90MEAN <-mean(cold90$n)
coldCI90 <- cold90SD/sqrt(100) 

herd_equal_Mean <- c(100-cold0MEAN,100-cold10MEAN,100-cold20MEAN,100-cold30MEAN ,100-cold40MEAN,
                       100-cold50MEAN,100-cold60MEAN,
                       100-cold70MEAN,100-cold80MEAN, 100-cold90MEAN)
herd_equal_SD <- c(cold0SD,cold10SD,cold20SD,cold30SD,cold40SD,cold50SD,cold60SD,cold70SD,cold80SD,cold90SD)
herd_equal_CI <- c(coldCI0,coldCI10,coldCI20,coldCI30,coldCI40,coldCI50,coldCI60,coldCI70,coldCI80,coldCI90)

Percentage_herd_immune <- c(0,10,20,30,40,50,60,70,80,90)
herd_equal_Mean[is.na(herd_equal_Mean)]= 100
herd_equal_SD[is.na(herd_equal_SD)]= 0
herd_equal_CI[is.na(herd_equal_CI)]= 0
herd_immunity_equal <- data.frame(herd_equal_Mean,herd_equal_SD ,herd_equal_CI)
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
########################################################################################################################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

endT = 30
deltaT = 1
time = seq(0,endT,by = deltaT)
temp = matrix(9,nrow = length(time))

beta <-matrix(0.0184,nrow = length(time))
sigma <- matrix(0.04199*temp+0.34936,nrow = length(time))
alpha <- matrix(-0.024049*temp+0.429996,nrow = length(time))
gamma <- matrix(0.058687*temp-0.453518,nrow = length(time))
tau <- matrix(0.01895*temp-0.04254,nrow = length(time))
delta <- matrix(-0.3959*temp+11.0789,nrow = length(time))
recoinf <- matrix(0.00184, nrow = length(time))
pie <- data.frame(recoinf, temp)
pie$recoinf <- ifelse(pie$temp >= 9.1, 0, pie$recoinf)
pie <- pie$recoinf
recoshed <- matrix(3.5, nrow = length(time))
f <- data.frame(recoshed, temp)
f$recoshed <- ifelse(f$temp >= 9.1, 0, f$recoshed)
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
### 0% of recovered in the population
#############
N = 100
E = 0
I = 1
C = N*0.0
P = 0
D = 0
S=N-C-I-E

iter = 100
iterations = 99

coldsim02 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                        "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                        "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim0 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim0 <- data.frame(coldsim0)
  
  coldsim0$iteration =0
  
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim01 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim01<-data.frame(coldsim01)
    coldsim01$iteration<-j
    coldsim0<-rbind(coldsim01, coldsim0)
  }
  coldsim0$iter <- i
  coldsim02 <- rbind(coldsim02, coldsim0)
}


#############
### 10% of recovered in the population
#############
N = 100
E = 0
I = 1
C = N*0.10
P = 0
D = 0
S=N-C-I-E

iter = 100

coldsim102 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim10 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim10 <- data.frame(coldsim10)
  
  coldsim10$iteration =0
  
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim101 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim101<-data.frame(coldsim101)
    coldsim101$iteration<-j
    coldsim10<-rbind(coldsim101, coldsim10)
  }
  coldsim10$iter <- i
  coldsim102 <- rbind(coldsim102, coldsim10)
}


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
iter = 100

coldsim202 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim20 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim20 <- data.frame(coldsim20)
  
  coldsim20$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim201 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim201<-data.frame(coldsim201)
    coldsim201$iteration<-j
    coldsim20<-rbind(coldsim201, coldsim20)
  }
  coldsim20$iter <- i
  coldsim202 <- rbind(coldsim202, coldsim20)
}

#############
### 30% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.30
P = 0
D = 0
S=N-C-I-E

coldsim302 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim30 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim30 <- data.frame(coldsim30)
  
  coldsim30$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim301 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim301<-data.frame(coldsim301)
    coldsim301$iteration<-j
    coldsim30<-rbind(coldsim301, coldsim30)
  }
  coldsim30$iter <- i
  coldsim302 <- rbind(coldsim302, coldsim30)
}


#############
### 40% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.40
P = 0
D = 0
S=N-C-I-E

coldsim402 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim40 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim40 <- data.frame(coldsim40)
  
  coldsim40$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim401 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim401<-data.frame(coldsim401)
    coldsim401$iteration<-j
    coldsim40<-rbind(coldsim401, coldsim40)
  }
  coldsim40$iter <- i
  coldsim402 <- rbind(coldsim402, coldsim40)
}



#############
### 50% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.50
P = 0
D = 0
S=N-C-I-E

coldsim502 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim50 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim50 <- data.frame(coldsim50)
  
  coldsim50$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim501 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim501<-data.frame(coldsim501)
    coldsim501$iteration<-j
    coldsim50<-rbind(coldsim501, coldsim50)
  }
  coldsim50$iter <- i
  coldsim502 <- rbind(coldsim502, coldsim50)
}


#############
### 60% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.60
P = 0
D = 0
S=N-C-I-E

coldsim602 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim60 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim60 <- data.frame(coldsim60)
  
  coldsim60$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim601 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim601<-data.frame(coldsim601)
    coldsim601$iteration<-j
    coldsim60<-rbind(coldsim601, coldsim60)
  }
  coldsim60$iter <- i
  coldsim602 <- rbind(coldsim602, coldsim60)
}

#############
### 70% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.70
P = 0
D = 0
S=N-C-I-E

coldsim702 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim70 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim70 <- data.frame(coldsim70)
  
  coldsim70$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim701 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim701<-data.frame(coldsim701)
    coldsim701$iteration<-j
    coldsim70<-rbind(coldsim701, coldsim70)
  }
  coldsim70$iter <- i
  coldsim702 <- rbind(coldsim702, coldsim70)
}


#############
### 80% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.80
P = 0
D = 0
S=N-C-I-E

coldsim802 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim80 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim80 <- data.frame(coldsim80)
  
  coldsim80$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim801 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim801<-data.frame(coldsim801)
    coldsim801$iteration<-j
    coldsim80<-rbind(coldsim801, coldsim80)
  }
  coldsim80$iter <- i
  coldsim802 <- rbind(coldsim802, coldsim80)
}



#############
### 90% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.90
P = 0
D = 0
S=N-C-I-E

coldsim902 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim90 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim90 <- data.frame(coldsim90)
  
  coldsim90$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim901 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim901<-data.frame(coldsim901)
    coldsim901$iteration<-j
    coldsim90<-rbind(coldsim901, coldsim90)
  }
  coldsim90$iter <- i
  coldsim902 <- rbind(coldsim902, coldsim90)
}

cold0<- coldsim02 %>% filter(coldsim02$time == 30)
cold0 <- cold0 %>%
  filter(S > 78.5) %>%
  count(iteration)

cold10<- coldsim102 %>% filter(coldsim102$time == 30)
cold10 <- cold10 %>%
  filter(S > 71) %>%
  count(iteration)

cold20<- coldsim202 %>% filter(coldsim202$time == 30)
cold20 <- cold20 %>%
  filter(S > 63) %>%
  count(iteration)

cold30<- coldsim302 %>% filter(coldsim302$time == 30)
cold30 <- cold30 %>%
  filter(S > 55) %>%
  count(iteration)

cold40<- coldsim402 %>% filter(coldsim402$time == 30)
cold40 <- cold40 %>%
  filter(S > 47) %>%
  count(iteration)

cold50<- coldsim502 %>% filter(coldsim502$time == 30)
cold50 <- cold50 %>%
  filter(S > 39) %>%
  count(iteration)

cold60<- coldsim602 %>% filter(coldsim602$time == 30)
cold60 <- cold60 %>%
  filter(S > 31) %>%
  count(iteration)

cold70<- coldsim702 %>% filter(coldsim702$time == 30)
cold70 <- cold70 %>%
  filter(S > 23) %>%
  count(iteration)

cold80<- coldsim802 %>% filter(coldsim802$time == 30)
cold80 <- cold80 %>%
  filter(S > 15) %>%
  count(iteration)

cold90<- coldsim902 %>% filter(coldsim902$time == 3)
cold90 <- cold90 %>%
  filter(S > 7) %>%
  count(iteration)


cold0SD <-sd(cold0$n)
cold0MEAN <-mean(cold0$n)
coldCI0 <- cold0SD/sqrt(100) 

cold10SD <-sd(cold10$n)
cold10MEAN <-mean(cold10$n)
coldCI10 <- cold10SD/sqrt(100)

cold20SD <-sd(cold20$n)
cold20MEAN <-mean(cold20$n)
coldCI20 <- cold20SD/sqrt(100) 

cold30SD <-sd(cold30$n)
cold30MEAN <-mean(cold30$n)
coldCI30 <- cold30SD/sqrt(100) 

cold40SD <-sd(cold40$n)
cold40MEAN <-mean(cold40$n)
coldCI40 <- cold40SD/sqrt(100) 

cold50SD <-sd(cold50$n)
cold50MEAN <-mean(cold50$n)
coldCI50 <- cold50SD/sqrt(100) 

cold60SD <-sd(cold60$n)
cold60MEAN <-mean(cold60$n)
coldCI60 <- cold60SD/sqrt(100) 

cold70SD <-sd(cold70$n)
cold70MEAN <-mean(cold70$n)
coldCI70 <- cold70SD/sqrt(100) 

cold80SD <-sd(cold80$n)
cold80MEAN <-mean(cold80$n)
coldCI80 <- cold80SD/sqrt(100) 

cold90SD <-sd(cold90$n)
cold90MEAN <-mean(cold90$n)
coldCI90 <- cold90SD/sqrt(100) 

herd_tenth_Mean <- c(100-cold0MEAN,100-cold10MEAN,100-cold20MEAN,100-cold30MEAN ,100-cold40MEAN,
                     100-cold50MEAN,100-cold60MEAN,
                     100-cold70MEAN,100-cold80MEAN, 100-cold90MEAN)
herd_tenth_SD <- c(cold0SD,cold10SD,cold20SD,cold30SD,cold40SD,cold50SD,cold60SD,cold70SD,cold80SD,cold90SD)
herd_tenth_CI <- c(coldCI0,coldCI10,coldCI20,coldCI30,coldCI40,coldCI50,coldCI60,coldCI70,coldCI80,coldCI90)

Percentage_herd_immune <- c(0,10,20,30,40,50,60,70,80,90)

herd_tenth_Mean[is.na(herd_tenth_Mean)]= 100
herd_tenth_SD[is.na(herd_tenth_SD)]= 0
herd_tenth_CI[is.na(herd_tenth_CI)]= 0
herd_immunity_tenth <- data.frame(herd_tenth_Mean,herd_tenth_SD ,herd_tenth_CI)


####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
########################################################################################################################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

endT = 30
deltaT = 1
time = seq(0,endT,by = deltaT)
temp = matrix(9,nrow = length(time))

beta <-matrix(0.0184,nrow = length(time))
sigma <- matrix(0.04199*temp+0.34936,nrow = length(time))
alpha <- matrix(-0.024049*temp+0.429996,nrow = length(time))
gamma <- matrix(0.058687*temp-0.453518,nrow = length(time))
tau <- matrix(0.01895*temp-0.04254,nrow = length(time))
delta <- matrix(-0.3959*temp+11.0789,nrow = length(time))
recoinf <- matrix(0.000184, nrow = length(time))
pie <- data.frame(recoinf, temp)
pie$recoinf <- ifelse(pie$temp >= 9.1, 0, pie$recoinf)
pie <- pie$recoinf
recoshed <- matrix(3.5, nrow = length(time))
f <- data.frame(recoshed, temp)
f$recoshed <- ifelse(f$temp >= 9.1, 0, f$recoshed)
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
### 0% of recovered in the population
#############
N = 100
E = 0
I = 1
C = N*0.0
P = 0
D = 0
S=N-C-I-E

iter = 100
iterations = 99

coldsim02 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                        "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                        "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim0 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim0 <- data.frame(coldsim0)
  
  coldsim0$iteration =0
  
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim01 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim01<-data.frame(coldsim01)
    coldsim01$iteration<-j
    coldsim0<-rbind(coldsim01, coldsim0)
  }
  coldsim0$iter <- i
  coldsim02 <- rbind(coldsim02, coldsim0)
}


#############
### 10% of recovered in the population
#############
N = 100
E = 0
I = 1
C = N*0.10
P = 0
D = 0
S=N-C-I-E

iter = 100

coldsim102 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim10 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim10 <- data.frame(coldsim10)
  
  coldsim10$iteration =0
  
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim101 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim101<-data.frame(coldsim101)
    coldsim101$iteration<-j
    coldsim10<-rbind(coldsim101, coldsim10)
  }
  coldsim10$iter <- i
  coldsim102 <- rbind(coldsim102, coldsim10)
}


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
iter = 100

coldsim202 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim20 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim20 <- data.frame(coldsim20)
  
  coldsim20$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim201 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim201<-data.frame(coldsim201)
    coldsim201$iteration<-j
    coldsim20<-rbind(coldsim201, coldsim20)
  }
  coldsim20$iter <- i
  coldsim202 <- rbind(coldsim202, coldsim20)
}

#############
### 30% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.30
P = 0
D = 0
S=N-C-I-E

coldsim302 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim30 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim30 <- data.frame(coldsim30)
  
  coldsim30$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim301 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim301<-data.frame(coldsim301)
    coldsim301$iteration<-j
    coldsim30<-rbind(coldsim301, coldsim30)
  }
  coldsim30$iter <- i
  coldsim302 <- rbind(coldsim302, coldsim30)
}


#############
### 40% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.40
P = 0
D = 0
S=N-C-I-E

coldsim402 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim40 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim40 <- data.frame(coldsim40)
  
  coldsim40$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim401 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim401<-data.frame(coldsim401)
    coldsim401$iteration<-j
    coldsim40<-rbind(coldsim401, coldsim40)
  }
  coldsim40$iter <- i
  coldsim402 <- rbind(coldsim402, coldsim40)
}



#############
### 50% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.50
P = 0
D = 0
S=N-C-I-E

coldsim502 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim50 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim50 <- data.frame(coldsim50)
  
  coldsim50$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim501 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim501<-data.frame(coldsim501)
    coldsim501$iteration<-j
    coldsim50<-rbind(coldsim501, coldsim50)
  }
  coldsim50$iter <- i
  coldsim502 <- rbind(coldsim502, coldsim50)
}


#############
### 60% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.60
P = 0
D = 0
S=N-C-I-E

coldsim602 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim60 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim60 <- data.frame(coldsim60)
  
  coldsim60$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim601 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim601<-data.frame(coldsim601)
    coldsim601$iteration<-j
    coldsim60<-rbind(coldsim601, coldsim60)
  }
  coldsim60$iter <- i
  coldsim602 <- rbind(coldsim602, coldsim60)
}

#############
### 70% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.70
P = 0
D = 0
S=N-C-I-E

coldsim702 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim70 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim70 <- data.frame(coldsim70)
  
  coldsim70$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim701 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim701<-data.frame(coldsim701)
    coldsim701$iteration<-j
    coldsim70<-rbind(coldsim701, coldsim70)
  }
  coldsim70$iter <- i
  coldsim702 <- rbind(coldsim702, coldsim70)
}


#############
### 80% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.80
P = 0
D = 0
S=N-C-I-E

coldsim802 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim80 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim80 <- data.frame(coldsim80)
  
  coldsim80$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim801 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim801<-data.frame(coldsim801)
    coldsim801$iteration<-j
    coldsim80<-rbind(coldsim801, coldsim80)
  }
  coldsim80$iter <- i
  coldsim802 <- rbind(coldsim802, coldsim80)
}



#############
### 90% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.90
P = 0
D = 0
S=N-C-I-E

coldsim902 <- data.frame("time" = integer(), "S" = integer(), "E" = integer(),
                         "I" = integer(), "C" = integer(), "P" = integer(),"D"= integer(),
                         "iterations" = integer())

for(i in 1:iter){
  
  inits = c(S =S, E =E, I = I, C = C , P = P, D = D)############################## Initial values
  paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
            param$delta, param$phi, param$tau, param$rho, param$pie)####################### Parameter values
  coldsim90 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim90 <- data.frame(coldsim90)
  
  coldsim90$iteration =0
  
  for(j in 1:iterations){
    inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
    coldsim901 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
    coldsim901<-data.frame(coldsim901)
    coldsim901$iteration<-j
    coldsim90<-rbind(coldsim901, coldsim90)
  }
  coldsim90$iter <- i
  coldsim902 <- rbind(coldsim902, coldsim90)
}

cold0<- coldsim02 %>% filter(coldsim02$time == 30)
cold0 <- cold0 %>%
  filter(S > 78.5) %>%
  count(iteration)

cold10<- coldsim102 %>% filter(coldsim102$time == 30)
cold10 <- cold10 %>%
  filter(S > 71) %>%
  count(iteration)

cold20<- coldsim202 %>% filter(coldsim202$time == 30)
cold20 <- cold20 %>%
  filter(S > 63) %>%
  count(iteration)

cold30<- coldsim302 %>% filter(coldsim302$time == 30)
cold30 <- cold30 %>%
  filter(S > 55) %>%
  count(iteration)

cold40<- coldsim402 %>% filter(coldsim402$time == 30)
cold40 <- cold40 %>%
  filter(S > 47) %>%
  count(iteration)

cold50<- coldsim502 %>% filter(coldsim502$time == 30)
cold50 <- cold50 %>%
  filter(S > 39) %>%
  count(iteration)

cold60<- coldsim602 %>% filter(coldsim602$time == 30)
cold60 <- cold60 %>%
  filter(S > 31) %>%
  count(iteration)

cold70<- coldsim702 %>% filter(coldsim702$time == 30)
cold70 <- cold70 %>%
  filter(S > 23) %>%
  count(iteration)

cold80<- coldsim802 %>% filter(coldsim802$time == 30)
cold80 <- cold80 %>%
  filter(S > 15) %>%
  count(iteration)

cold90<- coldsim902 %>% filter(coldsim902$time == 3)
cold90 <- cold90 %>%
  filter(S > 7) %>%
  count(iteration)


cold0SD <-sd(cold0$n)
cold0MEAN <-mean(cold0$n)
coldCI0 <- cold0SD/sqrt(100) 

cold10SD <-sd(cold10$n)
cold10MEAN <-mean(cold10$n)
coldCI10 <- cold10SD/sqrt(100)

cold20SD <-sd(cold20$n)
cold20MEAN <-mean(cold20$n)
coldCI20 <- cold20SD/sqrt(100) 

cold30SD <-sd(cold30$n)
cold30MEAN <-mean(cold30$n)
coldCI30 <- cold30SD/sqrt(100) 

cold40SD <-sd(cold40$n)
cold40MEAN <-mean(cold40$n)
coldCI40 <- cold40SD/sqrt(100) 

cold50SD <-sd(cold50$n)
cold50MEAN <-mean(cold50$n)
coldCI50 <- cold50SD/sqrt(100) 

cold60SD <-sd(cold60$n)
cold60MEAN <-mean(cold60$n)
coldCI60 <- cold60SD/sqrt(100) 

cold70SD <-sd(cold70$n)
cold70MEAN <-mean(cold70$n)
coldCI70 <- cold70SD/sqrt(100) 

cold80SD <-sd(cold80$n)
cold80MEAN <-mean(cold80$n)
coldCI80 <- cold80SD/sqrt(100) 

cold90SD <-sd(cold90$n)
cold90MEAN <-mean(cold90$n)
coldCI90 <- cold90SD/sqrt(100) 

herd_hundredth_Mean <- c(100-cold0MEAN,100-cold10MEAN,100-cold20MEAN,100-cold30MEAN ,100-cold40MEAN,
                     100-cold50MEAN,100-cold60MEAN,
                     100-cold70MEAN,100-cold80MEAN, 100-cold90MEAN)
herd_hundredth_SD <- c(cold0SD,cold10SD,cold20SD,cold30SD,cold40SD,cold50SD,cold60SD,cold70SD,cold80SD,cold90SD)
herd_hundredth_CI <- c(coldCI0,coldCI10,coldCI20,coldCI30,coldCI40,coldCI50,coldCI60,coldCI70,coldCI80,coldCI90)

Percentage_herd_immune <- c(0,10,20,30,40,50,60,70,80,90)

herd_hundredth_Mean[is.na(herd_hundredth_Mean)]= 100
herd_hundredth_SD[is.na(herd_hundredth_SD)]= 0
herd_hundredth_CI[is.na(herd_hundredth_CI)]= 0
herd_immunity_hundredth <- data.frame(herd_hundredth_Mean,herd_hundredth_SD ,herd_hundredth_CI)


herd_immunity <- data.frame(herd_immunity_equal,herd_immunity_tenth,herd_immunity_hundredth)
herd_immunity

linetype = 6



sensitivityplot <- ggplot(herd_immunity, aes(x = Percentage_herd_immune))+
  geom_ribbon(aes(ymin = herd_equal_Mean-(1.96 * herd_equal_CI), ymax = herd_equal_Mean+(1.96 * herd_equal_CI)), colour = "blue",size = 3)+
  geom_ribbon(aes(ymin = herd_tenth_Mean-(1.96 * herd_tenth_CI), ymax = herd_tenth_Mean+(1.96 * herd_tenth_CI)), colour = "darkblue",size = 3)+
  geom_ribbon(aes(ymin = herd_hundredth_Mean-(1.96 * herd_hundredth_CI), ymax = herd_hundredth_Mean+(1.96 * herd_hundredth_CI)), colour = "purple",size =3)+
  geom_line(aes(y = herd_equal_Mean, colour = "βC = βI"),linetype = 6, linewidth = 2)+
  geom_line(aes(y = herd_tenth_Mean, colour = "βC = βI*0.1"),linetype = 4, linewidth = 2)+
  geom_line(aes(y = herd_hundredth_Mean, colour = "βC = βI*0.01"),  linetype = 3, linewidth = 2)+
  ylab(label="Probability that an epizootic 
       occurs in 30 days")+
  xlab(label="")+  
  annotate("text", x = 10, y = 95, label = "A", size = 5)+
  scale_y_continuous(limits = c(0,100))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "white"))+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="plain"),
        legend.text = element_text(size = 20),
        legend.key.width = unit(1.7, 'cm'))+
  scale_colour_manual("",
                      breaks=c("βC = βI","βC = βI*0.1","βC = βI*0.01"),
                      values=c("blue4","black","purple4"))





png(filename = "Sensitivity_Analysis.png", height = 4500, width = 6400, res = 900 )
sensitivityplot

dev.off()

herd_immunity_equal
herd_immunity_hundredth
herd
herd_immunity_tenth

truc <- data.frame(coldsim90$time,coldsim90$I)
transformed_truc <- truc %>%
  pivot_longer(cols = c("coldsim90.I"))
meantruc <- ddply(transformed_truc, .variable = c("time","name"), summarise, values= mean(value))
plot(meantruc$time,meantruc$values)

