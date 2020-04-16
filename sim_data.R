library(ggplot)

setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Spring 2020/IBIO 831/Coral_symbiodinium_stats/")
n.areas <- 2
n.species <- 2
n.samples <- 20
n <- n.areas * n.species * n.samples

areas <-gl(n=n.areas,k=n/n.areas)
areas

species <-gl(n=n.species,k=n.samples,length=n)
species

depth <- runif(n=n, min=0, max=5)
depth_st <- (depth - mean(depth))/sd(depth)

hist(depth_st)

size <- rgamma(n = n,shape = 3, rate = 0.1)
size_st <- (size-mean(size))/sd(size)

hist(size_st)

Xmat <- model.matrix(~areas*species+depth_st+size_st)
Xmat

logit_funct <- function(x){
  return(log(x/(1-x)))
}

inv_logit_funct <- function(x){
  return(1/(1+exp(-x)))
}


beta.vec <- c(-0.5, -2.2, 3.5, -3, 2, -0.4)

lin.pred <- Xmat[,] %*%beta.vec  

res <- inv_logit_funct(lin.pred)

boxplot(res~species*areas)

ggplot()
plot(res~depth, col=areas,
     pch=as.character(factor(species, labels = c("A", "B"))))

plot(res~size, col=areas,
     pch=as.character(factor(species, labels = c("A", "B"))))

sim_data <- data.frame(clade_d_prop = res,
                       depth = depth,
                       size = size, 
                       species = species)

write.csv(sim_data, "simulated_data.csv")
