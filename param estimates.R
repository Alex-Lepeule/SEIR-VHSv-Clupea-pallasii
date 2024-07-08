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
##############################################################################################################
data <- read.csv("cumulative mortality.csv", header = TRUE, dec = ".", sep = ",")
endT = 25
deltaT = 1
time = seq(0,endT,by = deltaT)
data$time <- time
##############################################################################################################
##############################################################################################################
##############################################################################################################
########################################### Cold Temperature 8.8°C ###########################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################

# Deterministic model

seirmod = function(t, y, parms) { 
S = y[1]
E = y[2]
I = y[3]
C = y[4]
D = y[5]
with(as.list(parms), { 
  dS = -mu * S - beta * S * I
  dE = beta * S * I - (mu + sigma) * E
  dI = sigma * E - (mu + gamma+ alpha) * I
  dC = gamma * I - mu * C
  dD = mu*S + mu*E + (mu+ alpha)*I + mu*C
res = c(dS, dE, dI, dC, dD)
list(res)
})
}

# RSS calculation and fitting model

lfn=function(p){ times = seq(0, 25, by=1)
start = c(S=99, E=0, I=1, C = 0, D = 0)
paras=exp(c(mu=log(0.0008), beta = p[1],sigma=p[2], gamma=p[3], alpha=p[4]))
out1 = as.data.frame(ode(start, times=times,seirmod, paras))
n=length(data$mean.8C)
rss=sum((data$mean.8C-out1$D)^2)
return(log(rss)*(n/2)-n*(log(n)-log(2*pi)-1)/2)
}
plot(data$time,data$mean.8C)

# initial values

paras0 = log(c(0.05,0.2,0.07,0.24))
fit = optim(paras0, lfn, hessian = TRUE)

times = seq(0, 25, by=1)
paras = exp(c(mu = log(0.0008),
              beta = fit$par[1], sigma = fit$par[2],
              gamma = fit$par[3], alpha = fit$par[4]))
start = c(S=99, E=0, I=1, C = 0, D = 0)
out1 = as.data.frame(ode(start, times, seirmod, paras))
plot(out1$time, out1$D, xlab="Time", ylab="Cumulative mortality",
     type="l")
lines(data$time, data$mean.8C, col=2, type="l")
legend("bottomright", c("Paul's data",
                     "SEIR fit"), lty=c(1,1), col=c(2,1))

summary(lm(out1$D~data$mean.8C))
# MLEs:
round(exp(fit$par), 4)#approximation of parameters


##############################################################################################################
##############################################################################################################
##############################################################################################################
############################################ Ambient temperature 11.2°C ######################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################

# Deterministic model

seirmod = function(t, y, parms) { 
  S = y[1]
  E = y[2]
  I = y[3]
  C = y[4]
  D = y[5]
  with(as.list(parms), { 
    dS = -mu * S - beta * S * I
    dE = beta * S * I - (mu + sigma) * E
    dI = sigma * E - (mu + gamma+ alpha) * I
    dC = gamma * I - mu * C
    dD = mu*S + mu*E + (mu+ alpha)*I + mu*C
    res = c(dS, dE, dI, dC, dD)
    list(res)
  })
}

# RSS calculation and fitting model

lfn=function(p){ times = seq(0, 25, by=1)
start = c(S=99, E=0, I=1, C = 0, D = 0)
paras=exp(c(mu=log(0.0008), beta=p[1],sigma=p[2], gamma=p[3], alpha=p[4]))
out2 = as.data.frame(ode(start, times=times,seirmod, paras))
n=length(data$mean.11C)
rss=sum((data$mean.11C-out2$D)^2)
return(log(rss)*(n/2)-n*(log(n)-log(2*pi)-1)/2)
}
plot(data$time,data$mean.11C)

# initial values

paras0 = log(c(0.05, 0.62,0.15,0.11))
fit = optim(paras0, lfn, hessian = TRUE)

times = seq(0, 25, by=1)
paras = exp(c(mu = log(0.0008),
              beta = fit$par[1], sigma = fit$par[2],
              gamma = fit$par[3], alpha = fit$par[4]))
start = c(S=99, E=0, I=1, C = 0, D = 0)
out2 = as.data.frame(ode(start, times, seirmod, paras))
plot(out2$time, out2$D, xlab="Time", ylab="Prevalence",
     type="l")
lines(data$time, data$mean.11C, col=2, type="l")
legend("bottomright", c("Paul's data",
                        "SEIR fit"), lty=c(1,1), col=c(2,1))

summary(lm(out2$D~data$mean.11C))
# MLEs:
round(exp(fit$par), 4)#approximation of parameters

##############################################################################################################
##############################################################################################################
##############################################################################################################
############################################# Warm temperature 14.7°C ########################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################

# Deterministic model

seirmod = function(t, y, parms) { 
  S = y[1]
  E = y[2]
  I = y[3]
  C = y[4]
  D = y[5]
  with(as.list(parms), { 
    dS = -mu * S - beta * S * I
    dE = beta * S * I - (mu + sigma) * E
    dI = sigma * E - (mu + gamma+ alpha) * I
    dC = gamma * I - mu * C
    dD = mu*S + mu*E + (mu+ alpha)*I + mu*C
    res = c(dS, dE, dI, dC, dD)
    list(res)
  })
}

# RSS calculation and fitting model

lfn=function(p){ times = seq(0, 25, by=1)
start = c(S=99, E=0, I=1, C = 0, D = 0)
paras=exp(c(mu=log(0.0008), beta=p[1],sigma=p[2], gamma=p[3], alpha=p[4]))
out3 = as.data.frame(ode(start, times=times,seirmod, paras))
n=length(data$mean..14C)
rss=sum((data$mean..14C-out3$D)^2)
return(log(rss)*(n/2)-n*(log(n)-log(2*pi)-1)/2)
}
plot(data$time,data$mean.11C)

# initial values
paras0 = log(c(0.05,0.82,0.21,0.045))
fit = optim(paras0, lfn, hessian = TRUE)

times = seq(0, 25, by=1)
paras = exp(c(mu = log(0.0008),
              beta = fit$par[1], sigma = fit$par[2],
              gamma = fit$par[3], alpha = fit$par[4]))
start = c(S=99, E=0, I=1, C = 0, D = 0)
out3 = as.data.frame(ode(start, times, seirmod, paras))
plot(out3$time, out3$D, xlab="Time", ylab="Prevalence",
     type="l")
lines(data$time, data$mean..14C, col=2, type="l")
legend("bottomright", c("Paul's data",
                        "SEIR fit"), lty=c(1,1), col=c(2,1))

summary(lm(out3$D~data$mean..14C))
# MLEs:
round(exp(fit$par), 4)#approximation of parameters

plot(out1$time, out1$D, xlab="Time", ylab="Prevalence",
     type="l")
lines(data$time, data$mean.8C, col=2, type="l")
lines(meanDcold$time, meanDcold$values, col=2, type="l")
legend("bottomright", c("Paul's data",
                        "SEIR fit"), lty=c(1,1), col=c(2,1))

plot(out2$time, out2$D, xlab="Time", ylab="Prevalence",
     type="l")
lines(data$time, data$mean.11C, col=2, type="l")
lines(meanDAmbient$time, meanDAmbient$values, col=2, type="l")
legend("bottomright", c("SEIR fit",
                        "Paul's data",
                        "Stochastic"), lty=c(1,1), col=c(3,1))


plot(out3$time, out3$D, xlab="Time", ylab="Prevalence",
     type="l")
lines(data$time, data$mean..14C, col=2, type="l")
lines(meanDWarm$time, meanDWarm$values, col=2, type="l")
legend("bottomright", c("Paul's data",
                        "SEIR fit"), lty=c(1,1), col=c(2,1))

summary(lm(out3$D~data$mean..14C))
summary(lm(out2$D~data$mean.11C))
summary(lm(out1$D~data$mean.8C))

likecoldset <- data.frame(out1$D,data$time, data$mean.8C, data$lower.SD.8C,data$upper.SD.8C)
likecoldset
likecold <- ggplot(likecoldset,aes(x=data.time))+
  geom_line(aes(y=out1.D, colour = "Deterministic model results"), alpha = 1, size = 20)+
  geom_point(aes(y=data.mean.8C, colour = "Hershberger et al Data"), size = 30)+
  geom_errorbar(aes(ymin = data.lower.SD.8C, ymax = data.upper.SD.8C), size = 10, colour = "red")+
  ylab(label="Cumulative mortality (%)")+
  xlab(label="Days post exposure")+
  scale_y_continuous(limits = c(0,100))+
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  theme(axis.text=element_text(size=175),
        axis.title=element_text(size=175,face="bold"))+
  theme(legend.justification=c(1,0), legend.position=c(0.85,0.85))+
  theme(legend.title=element_text(size=0,face="bold"),
        legend.background = element_rect(fill='#FFFFFF',
                                         size=100,linetype="solid"),
        legend.text=element_text(size=175),
        legend.key=element_rect(colour="#FFFFFF",
                                fill='#FFFFFF',
                                size=0,
                                linetype="solid"))+
  scale_colour_manual("",
                      breaks=c("Deterministic model results","Hershberger et al Data"),
                      values=c("blue", "red"))
likecold


likeambset <- data.frame(out2$D,data$time, data$mean.11C, data$lower.SD.11C,data$upper.SD.11C)
likeambset
likeamb <- ggplot(likeambset,aes(x=data.time))+
  geom_line(aes(y=out2.D), colour = "blue", alpha = 1, size = 20)+
  geom_point(aes(y=data.mean.11C), size = 30, colour = "red")+
  geom_errorbar(aes(ymin = data.lower.SD.11C, ymax = data.upper.SD.11C), size = 10, colour = "red")+
  ylab(label="")+
  xlab(label="Days post exposure")+
  scale_y_continuous(limits = c(0,100))+
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  theme(axis.text=element_text(size=175),
        axis.title=element_text(size=175,face="bold"))
likeamb


likewarmset <- data.frame(out3$D,data$time, data$mean..14C, data$lower.SD.14C,data$upper.SD.14C)
likewarmset
likewarm <- ggplot(likewarmset,aes(x=data.time))+
  geom_line(aes(y=out3.D), colour = "blue", alpha = 1, size = 20)+
  geom_point(aes(y=data.mean..14C), size = 30, colour = "red")+
  geom_errorbar(aes(ymin = data.lower.SD.14C, ymax = data.upper.SD.14C), size = 10, colour = "red")+
  ylab(label="")+
  xlab(label="Days post exposure")+
  scale_y_continuous(limits = c(0,100))+
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 6),
        panel.background = element_rect(fill = "white"))+
  theme(axis.text=element_text(size=175),
        axis.title=element_text(size=175,face="bold"))
likewarm


png(filename = "estimates.png", height = 4500, width = 12800, res = 100 )

(likecold|likeamb|likewarm)

dev.off()
