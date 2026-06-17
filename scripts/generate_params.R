#install.packages('here')
library(here)

# mat_p_sim <- cbind(rep(c(5, 10, 20, 40, 80), c(1, 2, 4, 10, 20)),
#                    rep(c(100, 50, 25, 10, 5), c(1, 2, 4, 10, 20)))
mat_p_sim <- cbind(rep(c(5), c(5)),
                   rep(c(1), c(5)))
colnames(mat_p_sim) <- c('p', 'nsim')
params <- expand.grid(
  family = c("poisson"),
  signal =  NA,
  row_id = seq_len(nrow(mat_p_sim)),
  r = c(4, 8, 16, 32, 64),
  c_val =  NA_real_,
  sigma2 = NA,
  globalTest = c(FALSE),
  test = c('Both'),
  Swap = TRUE,
  alpha = log(9/4),
  beta = log(2),
  quantile = 0.35,
  Corr = 1,
  Eigen = c(1, 2),
  Method = "Method1",
  Method_para = 2,
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
          row.names = FALSE,
          quote = FALSE)

