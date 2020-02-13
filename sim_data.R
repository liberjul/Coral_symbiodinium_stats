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

Xmat = model.matrix(~areas*species+depth_st+size_st)
Xmat
