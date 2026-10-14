save_checkpoint <- function(results,
                            job_id,
                            p,
                            r,
                            Eigen,
                            family,
                            path = "results") {
  
  save_path <- file.path(
    path,
    paste0("Eigen", Eigen, "_", family, "_", "axis_method", axis_method)
  )
  
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  saveRDS(
    results,
    file = file.path(
      save_path,
      paste0("job_", job_id, "_p_", p, "_r_", r, ".rds")
    )
  )
}