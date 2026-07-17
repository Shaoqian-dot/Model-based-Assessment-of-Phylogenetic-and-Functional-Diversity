# When globalTest == TRUE, only poisson and gaussian families are allowed;
# When globalTest == FALSE and Swap == TRUE, only poisson and binomial families are allow; 
# When globalTest == FALSE and Swap == FALSE, only poisson and gaussian families are allow.
run_simulation <- function(
    signal_type,
    family_type,
    globalTest,
    test,
    Swap,
    p,
    q,
    r,
    c_val,
    sigma2,
    seed,
    tree,
    alpha,
    beta,
    quantile,
    Corr,
    Eigen,
    axis_method
){
  save_checkpoint <- function(results,
                              job_id,
                              p,
                              r,
                              Eigen,
                              family,
                              path = "results") {
    
    save_path <- file.path(
      path,
      paste0("Eigen", Eigen, "_", family)
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
  
  grid <- tidyr::expand_grid(
    p = p,
    r = r,
    family_type = family_type,
    Eigen = Eigen
  ) |>
    dplyr::mutate(
      sim = dplyr::row_number(),
      seed_i = seed + sim - 1
    )
  
  res <- purrr::pmap_dfr(
    grid,
    function(p, r, family_type, Eigen, sim, seed_i) {

      set.seed(seed_i)

      dir.create("timing_logs", showWarnings = FALSE)
      timing_file <- file.path(
        "timing_logs",
        sprintf("timing_seed_%d.csv", seed)
      )

      timing <- system.time({

        DM_phy_func <- getPhyloMatrix(tree = tree, m = p)
        VC_phy_func <- 1 - DM_phy_func / max(DM_phy_func)

        I_p <- diag(p)
        J <- I_p - 1 / p * matrix(1, p, p)
        S_J <- t(J) %*% VC_phy_func %*% J

        eig <- spectral_decomp(VC_phy_func = S_J)
        V_J <- eig$P

        cat(
          "Running:",
          signal_type,
          "| family =", family_type,
          "| p =", p,
          "| r =", r,
          "| Eigen =", Eigen,
          "\n"
        )

        if (Swap) {

          NOPS <- switch(
            as.character(p),
            "5" = 1,
            "10" = 2,
            "20" = 4,
            "40" = 8,
            "80" = 16,
            "160" = 32,
            stop("Unknown p")
          )

          dat <- getData(
            DM_phy_func = DM_phy_func,
            NOPS = NOPS,
            q = q,
            alpha = alpha,
            beta = beta,
            r = r,
            quantile = quantile,
            Corr = Corr,
            P = V_J[, 1:(p - 1)],
            Eigen = Eigen,
            Distribution = family_type,
            p = p
          )

          dat_long <- dat$yX
          dat_wide <- dat$abundance

        } else {

          dat <- generate_data(
            p = p,
            r = r,
            signal = signal_type,
            family = family_type,
            c_val = c_val,
            sigma2 = sigma2,
            tree = tree,
            V = V_J
          )

          dat_long <- dat$long
          dat_wide <- dat$wide[, -1]

        }

        model_res <- fit_models(
          dat = dat_long,
          family = family_type,
          globalTest = globalTest,
          tree = tree,
          S_J = S_J,
          p = p,
          q = q,
          test = test,
          axis_method = axis_method
        )

        # model_res <- model_res %>%
        #   tidyr::pivot_longer(
        #     cols = -c(axis_method, com, test),
        #     names_to = "model",
        #     values_to = "p_values"
        #   )
        # names(model_res)[names(model_res) %in% c("glmm_no_rr", "glmm_rr")] <-
        #   c("glmm_no_rr_pvalue", "glmm_rr_pvalue")
        # 
        # model_res <- model_res %>%
        #   tidyr::pivot_longer(
        #     cols = c(
        #       glmm_no_rr_pvalue,
        #       glmm_no_rr_status,
        #       glmm_rr_pvalue,
        #       glmm_rr_status
        #     ),
        #     names_to = c("model", ".value"),
        #     names_pattern = "(glmm_(?:no_rr|rr))_(pvalue|status)"
        #   ) %>%
        #   dplyr::mutate(method = "model_based")
        
        # names(model_res)[names(model_res) %in% c("glmm_no_rr", "glmm_rr")] <-
        #   c("glmm_no_rr_pvalue", "glmm_rr_pvalue")
        # model_res <- model_res %>%
        #   tidyr::pivot_longer(
        #     cols = c(
        #       glmm_no_rr_pvalue,
        #       glmm_no_rr_status,
        #       glmm_no_rr_warning_mags,
        #       glmm_no_rr_null_warning_msgs,
        #       glmm_rr_pvalue,
        #       glmm_rr_status,
        #       glmm_rr_warning_msgs,
        #       glmm_rr_null_warning_msgs
        #     ),
        #     names_to = c("model", ".value"),
        #     names_pattern = "(glmm_(?:no_rr|rr))_(pvalue|status|warning_type|null_warning_type)"
        #   ) %>%
        #   dplyr::rename(
        #     alt_warning_type = warning_type
        #   ) %>%
        #   dplyr::mutate(
        #     method = "model_based"
        #   )
        
        names(model_res)[names(model_res) %in% c("glmm_no_rr", "glmm_rr")] <-
          c("glmm_no_rr_pvalue", "glmm_rr_pvalue")
        
        model_res <- model_res %>%
          tidyr::pivot_longer(
            cols = c(
              glmm_no_rr_pvalue,
              glmm_no_rr_status,
              glmm_no_rr_warning_msgs,
              glmm_no_rr_null_warning_msgs,
              glmm_no_rr_convergence,
              glmm_no_rr_pdHess,
              glmm_rr_pvalue,
              glmm_rr_status,
              glmm_rr_warning_msgs,
              glmm_rr_null_warning_msgs,
              glmm_rr_convergence,
              glmm_rr_pdHess
            ),
            names_to = c("model", ".value"),
            names_pattern = "(glmm_(?:no_rr|rr))_(pvalue|status|warning_msgs|null_warning_msgs|convergence|pdHess)"
          ) %>%
          dplyr::rename(
            alt_warning_msgs = warning_msgs,
            null_warning_msgs = null_warning_msgs,
            alt_convergence = convergence,
            alt_pdHess = pdHess
          ) %>%
          dplyr::mutate(
            method = "model_based"
          )

        rao <- getRaosQ(
          abundance = dat_wide,
          DM_phy_func = DM_phy_func,
          use_randomization = 0,
          q = q
        )

        rao_rand <- getRaosQ(
          abundance = dat_wide,
          DM_phy_func = DM_phy_func,
          use_randomization = 1,
          q = q
        )

        # rao_res <- dplyr::bind_rows(
        #   tibble::tibble(
        #     com = 2:q,
        #     model = "RaosQ",
        #     pvalue = rao
        #   ),
        #   tibble::tibble(
        #     com = 2:q,
        #     model = "Randomized RaosQ",
        #     pvalue = rao_rand
        #   )
        # ) %>%
        #   dplyr::mutate(method = "raoQ")
        
        rao_res <- dplyr::bind_rows(
          tibble::tibble(
            com = 2:q,
            model = "RaosQ",
            pvalue = rao,
            status = "Success"
          ),
          tibble::tibble(
            com = 2:q,
            model = "Randomized RaosQ",
            pvalue = rao_rand,
            status = "Success"
          )
        ) %>%
          dplyr::mutate(method = "raoQ")

        result <- dplyr::bind_rows(
          model_res,
          rao_res
        ) %>%
          dplyr::mutate(
            seed = seed_i,
            signal = signal_type,
            family = family_type,
            p = p,
            r = r,
            c_val = c_val,
            globalTest = globalTest,
            Eigen = Eigen,
            .before = 1
          )

      })  # end system.time()

      timing_df <- tibble::tibble(
        seed = seed_i,
        signal = signal_type,
        family = family_type,
        p = p,
        r = r,
        Eigen = Eigen,
        elapsed_sec = timing["elapsed"]
      )

      write.table(
        timing_df,
        file = timing_file,
        sep = ",",
        row.names = FALSE,
        col.names = !file.exists(timing_file),
        append = TRUE
      )

      cat(
        sprintf(
          "Finished: family=%s p=%d r=%d Eigen=%d | %.2f hours\n",
          family_type,
          p,
          r,
          Eigen,
          timing["elapsed"] / 3600
        )
      )
      save_checkpoint(
        results = result,
        job_id = seed,
        p = p,
        r = r,
        Eigen = Eigen,
        family = family_type
      )
      return(result)
    }
  )
  res
}