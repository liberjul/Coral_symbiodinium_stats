library(ggplot)
library(ggpubr)
library(patchwork)

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

size <- exp(rnorm(n = n, mean = 1, sd = 0.4))
log_size <- log(size)
log_size_st <- (log_size-mean(log_size))/sd(log_size)

hist(log_size_st)

Xmat <- model.matrix(~areas*species+depth_st+log_size_st)
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


plot(res~depth, col=areas,
     pch=as.character(factor(species, labels = c("A", "B"))))

plot(res~size, col=areas,
     pch=as.character(factor(species, labels = c("A", "B"))))

sim_data <- data.frame(clade_d_prop = res,
                       depth = depth,
                       size = size, 
                       species = species,
                       area = areas)

write.csv(sim_data, "simulated_data.csv")

depth_plot <- ggplot(sim_data, aes(x = depth, y = res, color = areas, shape = species)) +
  geom_point() +
  labs(x = "Depth (m)",
       y = "Proportion of Symbiodium Clade D",
       color = "Collection Area",
       shape = "Coral Species") + 
  theme_pubr() +
  scale_shape_manual(values = 16:17,
                     breaks = 1:2,
                     labels = c("Pocillopora damicornis", "Montipora digitata"),
                     guide =
                       guide_legend(label.theme = element_text(size = 10, angle = 0, face = "italic"))) +
  theme(legend.position = "right")
depth_plot
size_plot <- ggplot(sim_data, aes(x = size, y = res, color = areas, shape = species)) +
  geom_point() +
  labs(x = bquote("Size"(~m^2)),
       y = "Proportion of Symbiodium Clade D",
       color = "Collection Area",
       shape = "Coral Species") + 
  theme_pubr() +
  scale_shape_manual(values = 16:17,
                     breaks = 1:2,
                     labels = c("Pocillopora damicornis", "Montipora digitata"),
                     guide =
                         guide_legend(label.theme = element_text(size = 10, angle = 0, face = "italic"))) +
  theme(legend.position = "right")
size_plot

g <- depth_plot + size_plot + plot_layout(guides = 'collect')
g
ggsave("Predictor_and_response_plots.png", g, width = 12, height = 8, units = "in")
