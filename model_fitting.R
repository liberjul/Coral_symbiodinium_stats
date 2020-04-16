setwd("~/Documents/Classes/SS 2020/IBIO 831/Projects/Coral_symbiodinium_stats")

library(bbmle)
library(ggplot2) 
library(GGally) 
library(tidyverse)
library(piecewiseSEM) 

d <- read.csv("simulated_data.csv", row.names = "X")
d
########Looking at the data###################
ggpairs(d)
hist(d$depth, xlab = "depth")
hist(d$size, xlab = "size")
hist(log(d$depth), xlab = "log(depth)")
###########Cleaning data################
d$depth_centered = d$depth - mean(d$depth)
d
