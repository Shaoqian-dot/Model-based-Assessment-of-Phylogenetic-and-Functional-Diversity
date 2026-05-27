#install.packages('here')
library(here)

mat_p_sim <- cbind(rep(c(2, 4, 8, 16), c(1, 2, 5, 10)),
                   rep(c(100, 50, 20, 10), c(1, 2, 5, 10)))
colnames(mat_p_sim) <- c('p', 'nsim')
params <- expand.grid(
  family = c("gaussian"),
  signal = c("v1", "vp"),
  row_id = seq_len(nrow(mat_p_sim)),
  r = 80,
  c_val = c(0.3, 0.6),
  sigma2 = 1,
  globalTest = TRUE,
  stringsAsFactors = FALSE
)

params <- cbind(
  params,
  mat_p_sim[params$row_id, ]
)

params$seed <- 10000 + seq_len(nrow(params))
params$row_id <- NULL

write.csv(params,
          here("config", "params.csv"),
          row.names = FALSE)
