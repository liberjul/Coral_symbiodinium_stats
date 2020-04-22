setwd("~/Documents/Classes/SS 2020/IBIO 831/Projects/Coral_symbiodinium_stats")

library(bbmle)
library(ggplot2) 
library(GGally) 
library(tidyverse)
library(piecewiseSEM) 

d <- read.csv("simulated_data.csv", row.names = "X")
d
########Looking at the data###################
d$species = as.factor(d$species)
d$area = as.factor(d$area)
ggpairs(d)
par(mfrow=c(3,2))
hist(d$depth, xlab = "depth")
hist(d$size, xlab = "size")
hist(log(d$depth), xlab = "log_depth")
hist(log(d$size), xlab = "log_size")
hist(log(d$size) - log(mean(d$size)), xlab = "log_cent_size")
hist(log(d$depth)-log(mean(d$depth)), xlab = "log_cent_depth")
dev.off()
###########Cleaning data################
#Transforming according to histograms
d$depth_log = log(d$depth)
d$size_cent = log(d$size)
d$size_cent = d$size_cent - mean(d$size_cent)
#####################################
d
par(mfrow=c(2,2))
plot(clade_d_prop ~ depth_log, data = d, col = species)
plot(clade_d_prop ~ size_cent, data = d, col = species)
plot(log(clade_d_prop) ~ depth_log, data = d, col = species)
plot(log(clade_d_prop) ~ size_cent, data = d, col = species)
dev.off()
#better represntation when we transform into log clade_d_prop
d$clade_d_prop_lg = log(d$clade_d_prop)
#################model selection########################
#additive
m0 <- lm(clade_d_prop_lg ~ 0, data = d)
m1 <- lm(clade_d_prop_lg ~ depth_log, data = d)
m2 <- lm(clade_d_prop_lg ~ size_cent, data = d)
m3 <- lm(clade_d_prop_lg ~ species, data = d)
m4 <- lm(clade_d_prop_lg ~ area, data = d)

m12 <- lm(clade_d_prop_lg ~ depth_log +size_cent, data = d)
m13 <- lm(clade_d_prop_lg ~ depth_log +species, data = d)
m14 <- lm(clade_d_prop_lg ~ depth_log +area, data = d)
m23 <- lm(clade_d_prop_lg ~ size_cent+species, data = d)
m24 <- lm(clade_d_prop_lg ~ size_cent+area, data = d)
m34 <- lm(clade_d_prop_lg ~ species+area, data = d)

m123 <- lm(clade_d_prop_lg ~ depth_log +size_cent + species , data = d)
m134 <- lm(clade_d_prop_lg ~ depth_log +size_cent + area , data = d)
m234 <- lm(clade_d_prop_lg ~ size_cent + species + area , data = d)

m1234 <- lm(clade_d_prop_lg ~ depth_log +size_cent + species + area , data = d)

# interative


m12_i <- lm(clade_d_prop_lg ~ depth_log *size_cent, data = d)
m13_i <- lm(clade_d_prop_lg ~ depth_log *species, data = d)
m14_i <- lm(clade_d_prop_lg ~ depth_log *area, data = d)
m23_i <- lm(clade_d_prop_lg ~ size_cent*species, data = d)
m24_i <- lm(clade_d_prop_lg ~ size_cent*area, data = d)
m34_i <- lm(clade_d_prop_lg ~ species*area, data = d)

m123_i <- lm(clade_d_prop_lg ~ depth_log *size_cent * species , data = d)
m134_i <- lm(clade_d_prop_lg ~ depth_log *size_cent * area , data = d)
m234_i <- lm(clade_d_prop_lg ~ size_cent * species * area , data = d)

m1234_i <- lm(clade_d_prop_lg ~ depth_log *size_cent * species * area , data = d)
#####just to see these make sense
m_fixed <- lmer(clade_d_prop_lg ~ depth_log +size_cent +species + (1|area) , data = d)
m_fixed_i <- lmer(clade_d_prop_lg ~ depth_log *size_cent *species + (1|area) , data = d)
m134_a <- lm(clade_d_prop_lg ~ depth_log +size_cent*area*species , data = d)
#########################################
tab = AICctab(m0, m1, m2, m3, m4, m12, m13, m14, m23, m24, m34,m123, m134,m234,m1234,
              m12_i, m13_i, m14_i, m23_i, m24_i, m34_i,m123_i, m134_i,m234_i,m1234_i,m_fixed,m_fixed_i,m134_a,
              base=TRUE, delta=TRUE, weights=TRUE, logLik=TRUE)
class(tab) = 'data.frame'
tab
#best model is m1234_i
##################looking at the model########################
summary(m1234_i)
rsquared(m1234_i)
#model explain 88.05% of the data
############################################
par(mfrow=c(2,2))
plot(residuals(m1234_i) ~ depth_log, data = d)
lines(lowess(residuals(m1234_i) ~ d$depth_log), col = 2)
plot(residuals(m1234_i) ~ size_cent, data = d)
lines(lowess(residuals(m1234_i) ~ d$size_cent), col = 2)
plot(residuals(m1234_i) ~ species, data = d)
plot(residuals(m1234_i) ~ area, data = d)
dev.off()

############################################
plot(resid(m1234_i) ~ predict(m1234_i))
lines(lowess(resid(m1234_i) ~ predict(m1234_i)), col=2)

########################################
summary(m1234_i)

newdata=data.frame(depth_log=seq(-3,3, length=1000), 
                   size_cent = seq(-2,2, length=1000),
                   species=as.factor(rep(c(1, 2), length=1000)),
                   area=as.factor(rep(c(1, 2), length=1000))
                   )

m1234_i_p = predict(m1234_i, newdata, se.fit=TRUE)

par(mfrow=c(2,2))
plot(clade_d_prop ~ depth, data = d, col = area)
lines(exp(m1234_i_p$fit) ~ exp(newdata$depth_log) , subset=newdata$area==1)
lines(exp(m1234_i_p$fit-m1234_i_p$se.fit*1.96) ~ exp(newdata$depth_log), subset=newdata$area==1, lwd=0.2)
lines(exp(m1234_i_p$fit+m1234_i_p$se.fit*1.96) ~ exp(newdata$depth_log), subset=newdata$area==1, lwd=0.2)
plot(clade_d_prop ~ depth, data = d, col = area)
lines(exp(m1234_i_p$fit) ~ exp(newdata$depth_log) , subset=newdata$area==2, col=2)
lines(exp(m1234_i_p$fit-m1234_i_p$se.fit*1.96) ~ exp(newdata$depth_log), subset=newdata$area==2, lwd=0.2, col = 2)
lines(exp(m1234_i_p$fit+m1234_i_p$se.fit*1.96) ~ exp(newdata$depth_log), subset=newdata$area==2, lwd=0.2, col = 2)




plot(clade_d_prop ~ size, data = d, col = area)
lines(exp(m1234_i_p$fit) ~ log(exp(newdata$size_cent+mean(d$size))) , subset=newdata$area==1)
lines(exp(m1234_i_p$fit-m1234_i_p$se.fit*1.96) ~ log(exp(newdata$size_cent+mean(d$size))), subset=newdata$area==1, lwd=0.2)
lines(exp(m1234_i_p$fit+m1234_i_p$se.fit*1.96) ~ log(exp(newdata$size_cent+mean(d$size))), subset=newdata$area==1, lwd=0.2)

plot(clade_d_prop ~ size, data = d, col = area)
lines(exp(m1234_i_p$fit) ~ log(exp(newdata$size_cent+mean(d$size))) , subset=newdata$area==2,col = 2)
lines(exp(m1234_i_p$fit-m1234_i_p$se.fit*1.96) ~ log(exp(newdata$size_cent+mean(d$size))), subset=newdata$area==2, lwd=0.2,col = 2)
lines(exp(m1234_i_p$fit+m1234_i_p$se.fit*1.96) ~ log(exp(newdata$size_cent+mean(d$size))), subset=newdata$area==2, lwd=0.2,col = 2)
dev.off()


