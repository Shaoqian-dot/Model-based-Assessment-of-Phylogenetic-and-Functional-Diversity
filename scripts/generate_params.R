#install.packages('here')
library(here)

############################################################
## Simulation settings
############################################################

nsim = 1000
#step <- length(p) * length(r) * length(families) * length(Eigen) 
step <- 100

params <- data.frame(
  seed = seq(
    from = 1000,
    by = step,
    length.out = nsim
  )
)

write.csv(params,
          here("config", "params.csv"),
          row.names = FALSE,
          quote = FALSE)

