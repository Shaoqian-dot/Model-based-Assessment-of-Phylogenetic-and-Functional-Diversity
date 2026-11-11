# Check if theta is factor loadings
## Extract factor loadings
report_obj <- Model_res_101$obj$env$report(Model_res_101$fit$parfull)
## Loadings for the first rr() term
Lambda <- report_obj$fact_load[[2]]

all.equal(c(Lambda), 
          unname(Model_res_101$fit$par[names(Model_res_101$fit$par) == 'theta']))



# Check if b is estimated random effects
## Extract estimated ramdom effects
source(system.file("misc", "extract_rr.R", package = "glmmTMB"))
rr_info <- extract_rr(Model_res_101)
RE <- rr_info$b %*% t(rr_info$fl)
b <- Model_res_101$fit$parfull[names(Model_res_101$fit$parfull) == 'b']
all.equal(c(RE),
          unname(b))
