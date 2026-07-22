files <- list.files(
  "output/PhyD_Pois_p80_r_TestBoth_MethodBoth_Eigen2",
  pattern = "\\.rds$",
  full.names = FALSE
)

job_ids <- unique(as.integer(
  sub("res_seed([0-9]+).*", "\\1", files)
))

expected <- 1000:1999

missing <- setdiff(expected, job_ids)

missing
length(missing)

missing + 1 - 1000
