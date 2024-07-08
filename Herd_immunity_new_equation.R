
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
           quote(ifelse(phi*S*C > S, S, phi*S*C)),
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
temp = matrix(8.8,nrow = length(time))

beta <-matrix(0.0184,nrow = length(time))
sigma <- matrix(0.04199*temp-0.34936,nrow = length(time))
alpha <- matrix(-0.024049*temp+0.429996,nrow = length(time))
gamma <- matrix(0.058687*temp-0.453518,nrow = length(time))
tau <- matrix(0.01895*temp-0.04254,nrow = length(time))
delta <- matrix(-0.3959*temp+11.0789,nrow = length(time))
recoshed <- matrix(0.0046, nrow = length(time))
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


SEICP = function(rateqs,eventmatrix, parameters, initialvals, deltaT, endT){
  time = seq(0,endT,by = deltaT)
  res <- data.frame(matrix(NA, ncol = length(initialvals)+ 1, nrow = length(time)))
  res[, 1] = time
  names(res) = c("time", names(inits))
  res[1, ] = c(0, inits)
  for(i in 1:(length(time)-1)){
    parameters <- c(beta = param$beta[i], sigma = param$sigma[i], gamma = param$gamma[i],
                    alpha = param$alpha[i], mu = param$mu[i], delta = param$delta[i],
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
coldsim0 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
coldsim0 <- data.frame(coldsim0)

coldsim0$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  coldsim01 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim01<-data.frame(coldsim01)
  coldsim01$iteration<-j
  coldsim0<-rbind(coldsim01, coldsim0)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
coldsim10 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
coldsim10 <- data.frame(coldsim10)

coldsim10$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  coldsim101 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim101<-data.frame(coldsim101)
  coldsim101$iteration<-j
  coldsim10<-rbind(coldsim101, coldsim10)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
coldsim20 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
coldsim20 <- data.frame(coldsim20)

coldsim20$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  coldsim201 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim201<-data.frame(coldsim201)
  coldsim201$iteration<-j
  coldsim20<-rbind(coldsim201, coldsim20)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
coldsim30 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
coldsim30 <- data.frame(coldsim30)

coldsim30$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  coldsim301 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim301<-data.frame(coldsim301)
  coldsim301$iteration<-j
  coldsim30<-rbind(coldsim301, coldsim30)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
coldsim40 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
coldsim40 <- data.frame(coldsim40)

coldsim40$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  coldsim401 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim401<-data.frame(coldsim401)
  coldsim401$iteration<-j
  coldsim40<-rbind(coldsim401, coldsim40)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
coldsim50 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
coldsim50 <- data.frame(coldsim50)

coldsim50$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  coldsim501 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim501<-data.frame(coldsim501)
  coldsim501$iteration<-j
  coldsim50<-rbind(coldsim501, coldsim50)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
coldsim60 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
coldsim60 <- data.frame(coldsim60)

coldsim60$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  coldsim601 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim601<-data.frame(coldsim601)
  coldsim601$iteration<-j
  coldsim60<-rbind(coldsim601, coldsim60)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
coldsim70 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
coldsim70 <- data.frame(coldsim70)

coldsim70$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  coldsim701 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim701<-data.frame(coldsim701)
  coldsim701$iteration<-j
  coldsim70<-rbind(coldsim701, coldsim70)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
coldsim80 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
coldsim80 <- data.frame(coldsim80)

coldsim80$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  coldsim801 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim801<-data.frame(coldsim801)
  coldsim801$iteration<-j
  coldsim80<-rbind(coldsim801, coldsim80)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
coldsim90 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
coldsim90 <- data.frame(coldsim90)

coldsim90$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  coldsim901 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim901<-data.frame(coldsim901)
  coldsim901$iteration<-j
  coldsim90<-rbind(coldsim901, coldsim90)
}

#############
### 95% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.95
P = 0
D = 0
S=N-C-I-E

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
coldsim95 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
coldsim95 <- data.frame(coldsim95)

coldsim95$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  coldsim951 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  coldsim951<-data.frame(coldsim951)
  coldsim951$iteration<-j
  coldsim95<-rbind(coldsim951, coldsim95)
}


cold0 <- subset(coldsim0, coldsim0$time == 365)
cold0 <- subset(cold0$iteration, cold0$S > 1)
cold0

cold10 <- subset(coldsim10, coldsim10$time == 365)
cold10 <- subset(cold10$iteration, cold10$S > 1)
cold10

cold20 <- subset(coldsim20, coldsim20$time == 365)
cold20 <- subset(cold20$iteration, cold20$S > 1)
cold20

cold30 <- subset(coldsim30, coldsim30$time == 365)
cold30 <- subset(cold30$iteration, cold30$S > 1)
cold30

cold40 <- subset(coldsim40, coldsim40$time == 365)
cold40 <- subset(cold40$iteration, cold40$S > 1)
cold40

cold50 <- subset(coldsim50, coldsim50$time == 365)
cold50 <- subset(cold50$iteration, cold50$S > 1)
cold50

cold60 <- subset(coldsim60, coldsim60$time == 365)
cold60 <- subset(cold60$iteration, cold60$S > 1)
cold60

cold70 <- subset(coldsim70, coldsim70$time == 365)
cold70 <- subset(cold70$iteration, cold70$S > 1)
cold70

cold80 <- subset(coldsim80, coldsim80$time == 365)
cold80 <- subset(cold80$iteration, cold80$S > 1)
cold80

cold90 <- subset(coldsim90, coldsim90$time == 365)
cold90 <- subset(cold90$iteration, cold90$S > 1)
cold90

cold95 <- subset(coldsim95, coldsim95$time == 365)
cold95 <- subset(cold95$iteration, cold95$S > 1)
cold95

num_epi_cold <- c(100-length(cold0),100-length(cold10),100-length(cold20),100-length(cold30),100-length(cold40),
                  100-length(cold50),100-length(cold60),
                  100-length(cold70),100-length(cold80),100-length(cold90),100-length(cold95))
Percentage_herd_immune <- c(0,10,20,30,40,50,60,70,80,90,95)

cold_herd_immunity <- data.frame(num_epi_cold,Percentage_herd_immune)

endT = 365
deltaT = 1
time = seq(0,endT,by = deltaT)
temp = matrix(11.2,nrow = length(time))

beta <-matrix(0.05,nrow = length(time))
sigma <- matrix(0.10159*temp-0.62839,nrow = length(time))
alpha <- matrix(-0.03413*temp+0.52245,nrow = length(time))
gamma <- matrix(0.019898*temp-0.086817,nrow = length(time))
tau <- matrix(0.01895*temp-0.04254,nrow = length(time))
delta <- matrix(-0.3959*temp+11.0789,nrow = length(time))
recoshed <- matrix(0.02, nrow = length(time))
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
ambientsim0 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
ambientsim0 <- data.frame(ambientsim0)

ambientsim0$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  ambientsim01 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  ambientsim01<-data.frame(ambientsim01)
  ambientsim01$iteration<-j
  ambientsim0<-rbind(ambientsim01, ambientsim0)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
ambientsim10 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
ambientsim10 <- data.frame(ambientsim10)

ambientsim10$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  ambientsim101 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  ambientsim101<-data.frame(ambientsim101)
  ambientsim101$iteration<-j
  ambientsim10<-rbind(ambientsim101, ambientsim10)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
ambientsim20 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
ambientsim20 <- data.frame(ambientsim20)

ambientsim20$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  ambientsim201 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  ambientsim201<-data.frame(ambientsim201)
  ambientsim201$iteration<-j
  ambientsim20<-rbind(ambientsim201, ambientsim20)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
ambientsim30 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
ambientsim30 <- data.frame(ambientsim30)

ambientsim30$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  ambientsim301 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  ambientsim301<-data.frame(ambientsim301)
  ambientsim301$iteration<-j
  ambientsim30<-rbind(ambientsim301, ambientsim30)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
ambientsim40 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
ambientsim40 <- data.frame(ambientsim40)

ambientsim40$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  ambientsim401 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  ambientsim401<-data.frame(ambientsim401)
  ambientsim401$iteration<-j
  ambientsim40<-rbind(ambientsim401, ambientsim40)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
ambientsim50 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
ambientsim50 <- data.frame(ambientsim50)

ambientsim50$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  ambientsim501 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  ambientsim501<-data.frame(ambientsim501)
  ambientsim501$iteration<-j
  ambientsim50<-rbind(ambientsim501, ambientsim50)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
ambientsim60 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
ambientsim60 <- data.frame(ambientsim60)

ambientsim60$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  ambientsim601 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  ambientsim601<-data.frame(ambientsim601)
  ambientsim601$iteration<-j
  ambientsim60<-rbind(ambientsim601, ambientsim60)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
ambientsim70 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
ambientsim70 <- data.frame(ambientsim70)

ambientsim70$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  ambientsim701 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  ambientsim701<-data.frame(ambientsim701)
  ambientsim701$iteration<-j
  ambientsim70<-rbind(ambientsim701, ambientsim70)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
ambientsim80 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
ambientsim80 <- data.frame(ambientsim80)

ambientsim80$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  ambientsim801 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  ambientsim801<-data.frame(ambientsim801)
  ambientsim801$iteration<-j
  ambientsim80<-rbind(ambientsim801, ambientsim80)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
ambientsim90 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
ambientsim90 <- data.frame(ambientsim90)

ambientsim90$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  ambientsim901 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  ambientsim901<-data.frame(ambientsim901)
  ambientsim901$iteration<-j
  ambientsim90<-rbind(ambientsim901, ambientsim90)
}


#############
### 95% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.95
P = 0
D = 0
S=N-C-I-E

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
ambientsim95 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
ambientsim95 <- data.frame(ambientsim95)

ambientsim95$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  ambientsim951 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  ambientsim951<-data.frame(ambientsim951)
  ambientsim951$iteration<-j
  ambientsim95<-rbind(ambientsim951, ambientsim95)
}

ambient0 <- subset(ambientsim0, ambientsim0$time == 365)
ambient0 <- subset(ambient0$iteration, ambient0$S > 1)
ambient0

ambient10 <- subset(ambientsim10, ambientsim10$time == 365)
ambient10 <- subset(ambient10$iteration, ambient10$S > 1)
ambient10

ambient20 <- subset(ambientsim20, ambientsim20$time == 365)
ambient20 <- subset(ambient20$iteration, ambient20$S > 1)
ambient20

ambient30 <- subset(ambientsim30, ambientsim30$time == 365)
ambient30 <- subset(ambient30$iteration, ambient30$S > 1)
ambient30

ambient40 <- subset(ambientsim40, ambientsim40$time == 365)
ambient40 <- subset(ambient40$iteration, ambient40$S > 1)
ambient40

ambient50 <- subset(ambientsim50, ambientsim50$time == 365)
ambient50 <- subset(ambient50$iteration, ambient50$S > 1)
ambient50

ambient60 <- subset(ambientsim60, ambientsim60$time == 365)
ambient60 <- subset(ambient60$iteration, ambient60$S > 1)
ambient60

ambient70 <- subset(ambientsim70, ambientsim70$time == 365)
ambient70 <- subset(ambient70$iteration, ambient70$S > 1)
ambient70

ambient80 <- subset(ambientsim80, ambientsim80$time == 365)
ambient80 <- subset(ambient80$iteration, ambient80$S > 1)
ambient80

ambient90 <- subset(ambientsim90, ambientsim90$time == 365)
ambient90 <- subset(ambient90$iteration, ambient90$S > 1)
ambient90

ambient95 <- subset(ambientsim95, ambientsim95$time == 365)
ambient95 <- subset(ambient95$iteration, ambient95$S > 1)
ambient95

num_epi_ambient <- c(100-length(ambient0),100-length(ambient10),100-length(ambient20),100-length(ambient30),100-length(ambient40),
                     100-length(ambient50),100-length(ambient60),
                     100-length(ambient70),100-length(ambient80),100-length(ambient90),100-length(ambient95))
Percentage_herd_immune <- c(0,10,20,30,40,50,60,70,80,90,95)

ambient_herd_immunity <- data.frame(num_epi_ambient,Percentage_herd_immune)


endT = 365
deltaT = 1
time = seq(0,endT,by = deltaT)
temp = matrix(14.7,nrow = length(time))

beta <-matrix(0.05,nrow = length(time))
sigma <- matrix(0.10159*temp-0.62839,nrow = length(time))
alpha <- matrix(-0.03413*temp+0.52245,nrow = length(time))
gamma <- matrix(0.019898*temp-0.086817,nrow = length(time))
tau <- matrix(0.01895*temp-0.04254,nrow = length(time))
delta <- matrix(-0.3959*temp+11.0789,nrow = length(time))
recoshed <- matrix(0.02, nrow = length(time))
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
warmsim0 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
warmsim0 <- data.frame(warmsim0)

warmsim0$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  warmsim01 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  warmsim01<-data.frame(warmsim01)
  warmsim01$iteration<-j
  warmsim0<-rbind(warmsim01, warmsim0)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
warmsim10 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
warmsim10 <- data.frame(warmsim10)

warmsim10$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  warmsim101 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  warmsim101<-data.frame(warmsim101)
  warmsim101$iteration<-j
  warmsim10<-rbind(warmsim101, warmsim10)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
warmsim20 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
warmsim20 <- data.frame(warmsim20)

warmsim20$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  warmsim201 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  warmsim201<-data.frame(warmsim201)
  warmsim201$iteration<-j
  warmsim20<-rbind(warmsim201, warmsim20)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
warmsim30 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
warmsim30 <- data.frame(warmsim30)

warmsim30$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  warmsim301 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  warmsim301<-data.frame(warmsim301)
  warmsim301$iteration<-j
  warmsim30<-rbind(warmsim301, warmsim30)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
warmsim40 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
warmsim40 <- data.frame(warmsim40)

warmsim40$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  warmsim401 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  warmsim401<-data.frame(warmsim401)
  warmsim401$iteration<-j
  warmsim40<-rbind(warmsim401, warmsim40)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
warmsim50 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
warmsim50 <- data.frame(warmsim50)

warmsim50$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  warmsim501 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  warmsim501<-data.frame(warmsim501)
  warmsim501$iteration<-j
  warmsim50<-rbind(warmsim501, warmsim50)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
warmsim60 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
warmsim60 <- data.frame(warmsim60)

warmsim60$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  warmsim601 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  warmsim601<-data.frame(warmsim601)
  warmsim601$iteration<-j
  warmsim60<-rbind(warmsim601, warmsim60)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
warmsim70 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
warmsim70 <- data.frame(warmsim70)

warmsim70$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  warmsim701 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  warmsim701<-data.frame(warmsim701)
  warmsim701$iteration<-j
  warmsim70<-rbind(warmsim701, warmsim70)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
warmsim80 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
warmsim80 <- data.frame(warmsim80)

warmsim80$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  warmsim801 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  warmsim801<-data.frame(warmsim801)
  warmsim801$iteration<-j
  warmsim80<-rbind(warmsim801, warmsim80)
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

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
warmsim90 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
warmsim90 <- data.frame(warmsim90)

warmsim90$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  warmsim901 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  warmsim901<-data.frame(warmsim901)
  warmsim901$iteration<-j
  warmsim90<-rbind(warmsim901, warmsim90)
}

#############
### 95% of recovered in the population
#############

N = 100
E = 0
I = 1
C = N*0.95
P = 0
D = 0
S=N-C-I-E

paras = c(param$beta, param$sigma, param$alpha, param$gamma, param$mu, 
          param$delta, param$phi, param$tau, param$rho)####################### Parameter values
inits = c(S =S, E =E, I = I, C = C , P = P, D = D)
warmsim95 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
warmsim95 <- data.frame(warmsim95)

warmsim95$iteration =0
iterations = 99

for(j in 1:iterations){
  inits = c(S = S, E =E, I = I, C = C, P = P, D = D)
  warmsim951 <- SEICP(rlist2, emat2, paras, inits,deltaT,endT)
  warmsim951<-data.frame(warmsim951)
  warmsim951$iteration<-j
  warmsim95<-rbind(warmsim951, warmsim95)
}


warm0 <- subset(warmsim0, warmsim0$time == 365)
warm0 <- subset(warm0$iteration, warm0$S > 1)
warm0

warm10 <- subset(warmsim10, warmsim10$time == 365)
warm10 <- subset(warm10$iteration, warm10$S > 1)
warm10

warm20 <- subset(warmsim20, warmsim20$time == 365)
warm20 <- subset(warm20$iteration, warm20$S > 1)
warm20

warm30 <- subset(warmsim30, warmsim30$time == 365)
warm30 <- subset(warm30$iteration, warm30$S > 1)
warm30

warm40 <- subset(warmsim40, warmsim40$time == 365)
warm40 <- subset(warm40$iteration, warm40$S > 1)
warm40

warm50 <- subset(warmsim50, warmsim50$time == 365)
warm50 <- subset(warm50$iteration, warm50$S > 1)
warm50

warm60 <- subset(warmsim60, warmsim60$time == 365)
warm60 <- subset(warm60$iteration, warm60$S > 1)
warm60

warm70 <- subset(warmsim70, warmsim70$time == 365)
warm70 <- subset(warm70$iteration, warm70$S > 1)
warm70

warm80 <- subset(warmsim80, warmsim80$time == 365)
warm80 <- subset(warm80$iteration, warm80$S > 1)
warm80

warm90 <- subset(warmsim90, warmsim90$time == 365)
warm90 <- subset(warm90$iteration, warm90$S > 1)
warm90

warm95 <- subset(warmsim95, warmsim95$time == 365)
warm95 <- subset(warm95$iteration, warm95$S > 1)
warm95

num_epi_warm <- c(100-length(warm0),100-length(warm10),100-length(warm20),100-length(warm30),100-length(warm40),
                  100-length(warm50),100-length(warm60),
                  100-length(warm70),100-length(warm80),100-length(warm90),100-length(warm95))
Percentage_herd_immune <- c(0,10,20,30,40,50,60,70,80,90,95)

warm_herd_immunity <- data.frame(num_epi_warm,Percentage_herd_immune)

herd_immunity <- data.frame(cold_herd_immunity,ambient_herd_immunity,warm_herd_immunity)

warmherdPlot <- ggplot(warm_herd_immunity,aes(x= Percentage_herd_immune))+
  geom_line(aes(y=num_epi_warm),size = 1, colour = "red")+
  ylab(label="Probability of an epizootic")+
  scale_y_continuous(labels = percent_format(scale = 1), limits = c(0,100))+
  scale_x_continuous(labels = percent_format(scale = 1), limits = c(0,100))+
  xlab(label="Percentage of carriers in the population")+
  theme_classic()+
  theme(axis.text=element_text(size=40),
        axis.title=element_text(size=40,face="bold"))
warmherdPlot

herdPlot <- ggplot(herd_immunity,aes(x= Percentage_herd_immune))+
  geom_line(aes(y=num_epi_cold),size = 3, colour = "blue4", linetype = "twodash")+
  geom_line(aes(y=num_epi_ambient),size = 3, colour = "orange", linetype = "dotdash")+
  geom_line(aes(y=num_epi_warm),size = 3, colour = "red", linetype = "dotted")+
  ylab(label="Probability of an epizootic")+
  scale_y_continuous(labels = percent_format(scale = 1), limits = c(0,100))+
  scale_x_continuous(labels = percent_format(scale = 1), limits = c(0,100))+
  xlab(label="Percentage of carriers in the population")+
  theme_classic()+
  theme(axis.text=element_text(size=40),
        axis.title=element_text(size=40,face="bold"))

png(filename = "Herd_Immunity.png", height = 900, width = 1280, res = 100 )
herdPlot

dev.off()
herdPlot
