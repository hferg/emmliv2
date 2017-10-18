###########################################################################################
#' identifyCandidates
#'
#' Identifies candidates for merging in simplifyEMMLi. Candidates are selected if their
#' between-module rho is larger than either of their within-module rhos by a factor of
#' 2*sd of the combined within-module correlations of the landmarks in the modules when
#' pooled.
#' @name simplifyEMMLi
#' @param fitted_emmli Fitted EMMLi model output.
#' @param corr_matrix The correlation matrix that EMMLi was fitted to.
#' @param models The original models that EMMLi tested.
#' @keywords internal

identifyCandidates <- function(fitted_emmli, corr_matrix, models) {
  fit_mods <- fitted_emmli$results
  best_mod <- names(fitted_emmli$rho)
  start_mod <- data.frame(data.point = models$Data.point,
                          original_best = models[ , grep(strsplit(best_mod, "\\.")[[1]][1], colnames(models))])
  rhos <- fitted_emmli$rho[[1]]

  # Seperate between from within.
  within_rhos <- rhos[ , grep("^Module", colnames(rhos))]
  between_rhos <- rhos[ , grep("to", colnames(rhos))]

  all_corrs <- getCorrs(fitted_emmli, models, corr_matrix)[[1]]
  between_rhos <- between_rhos[,order(between_rhos[2,], decreasing = TRUE)]

  # Select pairs that have a between rho greater than either of the withins by
  # a factor of 2*SD of the combined within-correlations of both modules pooled.

  pairs <- colnames(between_rhos)
  candidates <- list()
  for (i in seq_along(pairs)) {
    modules <- strsplit(pairs[1], " to ")[[1]]
    mod_1_rho <- rhos["MaxL_p" , grep(paste("Module", modules[1]), colnames(rhos))]
    mod_2_rho <- rhos["MaxL_p" , grep(paste("Module", modules[2]), colnames(rhos))]
    between_rho <- rhos["MaxL_p", grep(paste(modules[1], "to", modules[2]), colnames(rhos))]
    combined_mods <- c(
      all_corrs[[paste("Module", modules[1])]],
      all_corrs[[paste("Module", modules[2])]]
    )
    sd_comb <- sd(combined_mods)
    if (between_rho - mod_1_rho > 2 * sd_comb) {
      candidates[[length(candidates) + 1]] <- pairs[i]
    } else if (between_rho - mod_2_rho > 2 * sd_comb) {
      candidates[[length(candidates) + 1]] <- pairs[i]
    }
  }

  candidates <- lapply(candidates, function(x) as.numeric(strsplit(x, " to ")[[1]]))
  # return candidates, and the starting model (the best model from the fitted EMMLi)
  return(list(candidates = candidates, start_mod = start_mod))
}

###########################################################################################
#' simplifyEMMLi
#'
#' Attempts to simplify a fitted EMMLi model, looking for a simpler model that fits the
#' data better by merging modules. Pairs of modules are identified as candidates for merging
#' if their between-module rho is higher than either of the within-module rhos by more than
#' 2 * SD of the modules combined. This repeats until the model does not change or improve.
#' Alternatively, pairs of modules can be offered as candidates for merging - in this
#' instance there is not exploration and just the suggested pairs are tested. Pairs are
#' offered as a list, where each element is a vector of two module numbers.
#' @name simplifyEMMLi
#' @param fitted_emmli Fitted EMMLi model output.
#' @param corr_matrix The correlation matrix that EMMLi was fitted to.
#' @param models The original models that EMMLi tested.
#' @param candidates A list of pairs of modules to test merging (each element is a vector of
#' length 2 with module numbers).
#' @param N_sample The sample size for the original EMMLi fit.
#' @param correction if "normal" uses the normal EMMLi calculation for K, if "new" then uses
#' the experimental adjustment to AICc (adding nmodules - 1 to K). Not reccommended.
#' @export

simplifyEMMLi <- function(fitted_emmli, corr_matrix, models, candidates = NULL,
                          correction = "normal", N_sample) {

  if (is.null(candidates)) {
    x <- identifyCandidates(fitted_emmli, corr_matrix, models)
    candidates <- x$candidates
    start_mod <- start_mod
    if (length(candidates) == 0) {
      stop("No candidate pairs found.")
    }
  } else {
    # Use user-supplied candidates and identify starting model.
    candidates <- candidates
    best_mod <- names(fitted_emmli$rho)
    start_mod <- data.frame(data.point = models$Data.point,
                            original_best = models[ , grep(strsplit(best_mod, "\\.")[[1]][1], colnames(models))])
  }

  repeat {
    if (length(candidates) == 0) {
      break
    }

    new_emms <- list()
    # test all candidate pairs.
    for (i in seq_along(candidates)) {
      test_mods <- start_mod
      # Combine the two modules in a new_model
      test_mods$new_model <- test_mods$original_best
      test_mods$new_model[test_mods$new_model == candidates[[i]][1]] <- candidates[[i]][2]
      # fit EMMLi
      new_emms[[i]] <- EMMLi(corr = corr_matrix, mod = test_mods, N_sample = N_sample, correction = correction)

      # Now if the best model is "original" go to the next i, otherwise break this loop.
      # if (!grepl("original_best", rownames(emm$results)[emm$results[ , "dAICc"] == 0])) {
      #   print("Better model found.")
      #   models <- test_mods
      #   new_emm <- emm
      #   break
      # }

    }
    best_liks <- rep(NA, length(candidates))
    for (i in seq_along(new_emms)) {
      # If the original model isn't the best record the likelihood of it.
      if (!grepl("original_best",
                 rownames(new_emms[[1]]$results)[new_emms[[1]]$results[ , "dAICc"] == 0])) {
        best_liks[i] <- new_emms[[i]]$results[new_emms[[i]]$results[ , "dAICc"] == 0, "MaxL"]
      }
    }

    # If best_liks is all NA then none of the merges helped, and break the loop.
    # Else take the one with the best likelihood and calculate new canidate pairs.
    # If there are no candidates the loop breaks, if there are, it repeats.
    if (all(is.na(best_liks))) {
      break
    } else {
      new_emm <- new_emms[which.max(best_liks)]
      x <- identifyCandidates(fitted_emmli = new_emm, models = test_mods,
                              corr = corr_matrix)
      candidates <- x$candidates
      start_mod <- x$start_mod
    }
  }
  # If new_emm exists it will be the simplest version of the model, so return it
  # and the test_mods. This will be the most recent test mods.
  # If it doesn't say so, and break.
  if (exists("new_emm")) {
    return (list(emmli = new_emm, model = test_mods))
  } else {
    stop("No simplification possible.")
  }
}
