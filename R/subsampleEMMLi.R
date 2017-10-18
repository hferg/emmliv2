
###########################################################################################
#' subsampleLandmarks
#'
#' This function randomly removes landmarks from a dataset, and adjusts the corresponding
#' models object to match the newly subsampled landmark data. Internal, called by
#' subsampleEMMLi.
#' @name subsampleLandmarks
#' @param landmarks A 2D dataframe of landmarks
#' @param fraction The decimal fraction to subsample down to. e.g. 0.2 will return 20% of the original landmarks.
#' @param models The models for the landmarks
#' @param min_landmarks The minimum number of landmarks for a module in any model. When subsampling goes below this threshold landmarks are randomly seleclted from the original data and added back in in order to make up to min_landmarks. This means that sometimes the exact fraction desired is not returned.
#' @export
#' @keywords internal

subsampleLandmarks <- function(landmarks, fraction, models, min_landmark) {
  sampleSpecimen <- function(lms, kps) {
    return(lms[kps, ])
  }

  x <- lapply(models[2:length(models)], table)
  x <- sapply(x, length)
  kps <- sort(sample(nrow(landmarks[,,1]), round(nrow(landmarks[,,1]) * fraction, 0)))
  m <- models[ , 2:ncol(models)]
  for (i in 1:ncol(m)) {
    mod_original <- unique(m[ , i])
    ss_mod_count <- table(m[kps, i])

    if (length(ss_mod_count) != length(mod_original)) {
      missing <- mod_original[!mod_original %in% names(ss_mod_count)]
      for (j in missing) {
        kps <- c(kps, sample(which(m[ , i] == j), min_landmark))
      }
      kps <- sort(unique(kps))
    }

    if (any(ss_mod_count < min_landmark)) {
      short <- names(ss_mod_count)[ss_mod_count < min_landmark]
      for (j in short) {
        kps <- c(kps, sample(which(m[ , i] == j), min_landmark))
      }
      kps <- sort(unique(kps))
    }
  }

  res <- array(dim = c(length(kps), 3, dim(landmarks)[3]))
  for (i in 1:dim(landmarks)[3]) {
    res[,,i] <- sampleSpecimen(landmarks[,,i], kps)
  }

  new_models <- models[kps, ]
  return(list(landmarks = res, models = new_models,
              true_subsample = round(dim(res)[1] / dim(landmarks)[1], 2)))
}


###########################################################################################
#' subSampleEMMLi
#' This takes some landmarks, and then will subsample them down to one or more subsampling
#' fractions, calculate the correaltion matrix (using dotcorr only) and then fit EMMLi,
#' returning the result of all the subsampling (if multiple fractions). If a single fraction
#' the nsim argument is required with a number of times to repeat subsampling at that level.
#' @name subSampleEmmli
#' @param landmarks Landmarks dataframe
#' @param fractions Either a single subsampling fraction (in which case nsim is required) or a vector of fractions.
#' @param models The models to test.
#' @param min_landmark The minimum number of landmarks to subsample to
#' @param aic_cut This is the threshold of dAICc below which two models are considered to be not different. When this occurs multiple models can be returned as the best.
#' @param return_corr Logical - if TRUE then the full correlations of within and between modules are returned for the best models.
#' @param nsim If a single subsampling fraction, this is the number of times that subsampling fraction is repeated.
#' @export

subSampleEMMLi <- function(landmarks, fractions, models, min_landmark, aic_cut = 7,
                           return_corr = TRUE, nsim = NULL) {
  # Much like before, except this one takes actual data and actual models rather than
  # simulation parameters.
  require(paleomorph)
  require(tibble)
  #require(EMMLi)

  if (!is.null(nsim)) {
    if (length(fractions) > 1) {
      stop("Only one subsampling fraction is allowed when multiple subsampling simulations are run.")
    }
    fractions <- rep(fractions, nsim)
  }

  fitEmmli <- function(landmarks, models, fraction, min_landmark, aic_cut) {

    dat <- subsampleLandmarks(landmarks = landmarks, fraction = fraction, models = models,
                              min_landmark = min_landmark)
    c <- dotcorr(dat$landmarks)

    emm <- EMMLi(corr = c, mod = dat$models, N_sample = dim(landmarks)[3], all_rhos = TRUE)
    emmli <- as.data.frame(emm$results)
    best_models <- emmli[emmli$dAICc <= aic_cut, ]
    best_names <- rownames(best_models)
    best_models <- do.call(rbind, best_models)
    colnames(best_models) <- best_names
    names(emm$all_rhos) <- sapply(names(emm$all_rhos), function(x) strsplit(x, "\\$")[[1]][2])
    names(emm$all_rhos) <- unlist(strsplit(names(emm$all_rhos), "1$"))
    rhos <- emm$all_rhos[names(emm$all_rhos) %in% best_names]
    all_data <- list(landmarks = dat$landmarks, models = dat$models)

    if (return_corr) {
      best_corr <- getCorrs(emm, dat$models, c)
      res <- list(best_models = best_models, rhos_best = rhos, best_corr = best_corr,
                  all_data = all_data, true_subsample = dat$true_subsample)
    } else {
      res <- list(best_models = best_models, rhos_best = rhos, all_data = all_data,
                  true_subsample = dat$true_subsample)
    }
    # print(ncol(best_models))
    # if (ncol(best_models) > 1) {
    #   print("multiple best")
    # }
    return(res)
  }

  res <- lapply(fractions, function(x) fitEmmli(fraction = x, landmarks = landmarks,
                                                models = models, min_landmark = min_landmark, aic_cut = aic_cut))

  names(res) <- fractions

  return(res)
}
