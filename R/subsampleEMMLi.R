###########################################################################################
#' getCorrs
#'
#' This function takes a fitted EMMLi model and the models used to fit it, and then returns
#' all the correaltions for the best fitting model. This is recycled from the main EMMLi
#' function - Possible that this will change with the refactoring of EMMLi.
#' @name getCorrs
#' @param emm A fitted EMMLi model (output of EMMLi).
#' @param models A data frame defining the models. The first column should contain the landmark names
#' as factor or character. Subsequent columns should define which landmarks are contained within each
#' module with integers, factors or characters. If a landmark should be ignored for a specific model
#' (i.e., it is unintegrated in any module), the element should be NA.
#' @param corr The original correlation matrix used in the EMMLi analysis that generated
#' the input.
#' @return The correlations within and between modules of the best model.
#' @export

getCorrs <- function(emm, models, corr) {
  best_mod <- rownames(emm$results)[which(emm$res[ , "dAICc"] == 0)]
  best_mod <- strsplit(best_mod, " ")[[1]][1]
  best_mod <- paste(head(strsplit(best_mod, "\\.")[[1]], n = -2), collapse = ".")

  symmet = corr
  symmet[upper.tri(symmet)] = t(symmet)[upper.tri(symmet)]

  model <- paste0("models$", best_mod)
  lms <- array(eval(parse(text = model)))
  modNF <- stats::na.omit(cbind(1:nrow(models), lms))
  w <- unique(modNF[, 2])

  all_modules <- list()
  modules <- list()
  btw_mod = list()
  betweenModules = list()
  withinModules = list()

  for(i in seq(length(w))){
    # identify landmarks within a class
    fg <- modNF[modNF[, 2] == w[i], ]

    # coefficients between identified landmarks.
    l <- corr[fg[, 1], fg[, 1]]
    modules[[i]] <- (as.array(l[!is.na(l)]))
  }
  names(modules) <- paste("Module", w)

  if (length(w) > 1) {
    # make combinations of modules for between module.
    cb <- utils::combn(w, 2)
    for (i in seq(dim(cb)[2])){
      fg1 <- modNF[modNF[, 2] == cb[1, i], ]
      fg2 <- modNF[modNF[, 2] == cb[2, i], ]

      # setdiff(A,B) - present in A but not in B
      btw <- symmet[as.integer(setdiff(fg2[, 1], fg1[, 1])), as.integer(setdiff(fg1[, 1], fg2[, 1]))]
      btw_mod[[i]] <- btw[!is.na(btw)]
    }
    names(btw_mod) <- paste(cb[1, ], "to", cb[2, ])
    betweenModules['betweenModules'] = list(as.vector(rle(unlist(btw_mod))$values))
    withinModules['withinModules'] = list(as.vector(rle(unlist(modules))$values))
  }

  all_modules[[model]] = c(modules, btw_mod, withinModules, betweenModules)

  return(all_modules)
}


###########################################################################################
#' subsampleLandmarks
#'
#' This function randomly removes landmarks from a dataset, and adjusts the corresponding
#' models object to match the newly subsampled landmark data. Internal, called by
#' subsampleEMMLi.
#' @name subsampleLandmarks
#' @param landmarks A 2D array of xyz landmarks to subsample from. Will be turned into a correlation
#' matrix using \link[paleopmorph]{dotcorr} for EMMLi analysis after subsampling.
#' @param fraction The decimal fraction to subsample down to. e.g. 0.2 will return 20% of the
#' original landmarks.
#' @param models A data frame defining the models. The first column should contain the landmark names
#' as factor or character. Subsequent columns should define which landmarks are contained within each
#' module with integers, factors or characters. If a landmark should be ignored for a specific model
#' (i.e., it is unintegrated in any module), the element should be NA.
#' @param min_landmark The minimum number of landmarks to subsample to. This ensures that a module isn't
#' totally removed during random subsampling. When subsampling causes a landmark to be removed or to be
#' subsampled below this threshold landmarks are drawn from the original module at random and added back in.
#' This means that sometimes (especially with low subsampling fractions and/or the presence of small
#' modules in a model) the actual subsampling level is higher than the requested subsampling. In these
#' cases a warning is printed to the screen.
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

  true_subsample = round(dim(res)[1] / dim(landmarks)[1], 2)
  if (true_subsample > fraction) {
    warning(paste("Actual subsample is ", true_subsample))
  }

  new_models <- models[kps, ]
  return(list(landmarks = res, models = new_models, true_subsample = true_subsample))
}


###########################################################################################
#' subsampleEMMLi
#'
#' Analyse random subsamples of a dataset repeatedly using EMMLi. This function subsamples
#' a 2D array of landmarks to a given fraction of it's original size and then fits EMMLi
#' to the smaller dataset, returning the results. This can be either done repeatedly for
#' a single subsampling fraction, or for a range of subsampling fractions (e.g. to investigate
#' the effects of increasing subsampling). In all cases EMMLi is fit to a correlation matrix
#' calculated from the subsampled landmarks using \link[paleomorph]{dotcorr}.
#'
#' If a single fraction is provided then an nrep argument is also required, and the function
#' will subsample the data to the given fraction nrep times. Alternatively, if a range of
#' fractions is given (in a vector) the nrep argument is not required, and the function will
#' subsample the data and fit EMMLi once for each fraction in the given fractions vector.
#' @name subsampleEMMLi
#' @param landmarks A 2D array of xyz landmarks to subsample from. Will be turned into a correlation
#' matrix using \link[paleopmorph]{dotcorr} for EMMLi analysis after subsampling.
#' @param fractions Either a single subsampling fraction (in which case nrep is required) or a
#' vector of fractions. Specified in decimal format, i.e. a fraction of 0.2 will subsample down
#' to 20\% of the original number of landmarks.
#' @param models A data frame defining the models. The first column should contain the landmark names
#' as factor or character. Subsequent columns should define which landmarks are contained within each
#' module with integers, factors or characters. If a landmark should be ignored for a specific model
#' (i.e., it is unintegrated in any module), the element should be NA.
#' @param min_landmark The minimum number of landmarks to subsample to. This ensures that a module isn't
#' totally removed during random subsampling. When subsampling causes a landmark to be removed or to be
#' subsampled below this threshold landmarks are drawn from the original module at random and added back in.
#' This means that sometimes (especially with low subsampling fractions and/or the presence of small
#' modules in a model) the actual subsampling level is higher than the requested subsampling. In these
#' cases a warning is printed to the screen.
#' @param aic_cut This is the threshold of dAICc below which two models are considered to be not different.
#' When this occurs multiple models are be returned as the best. Defaults to 0 (i.e., only the best
#' fitting model is returned, regardless of how close other models are.)
#' @param return_corr Logical - if TRUE then the full correlations of within and between modules are
#' returned for the single best model in addition to the EMMLi results. Defaults to FALSE.
#' @param nrep If a single subsampling fraction, this is the number of times that subsampling fraction
#' is repeated.
#' @return A list of n elements, where n is either the number of subsampling fractions, or nrep. Each
#' element contains the output of an EMMLi analysis on the subasampled data, consisting of four or five
#' elements:
#'
#'   - the best model(s) description
#'
#'   - the rho list(s) for the best model(s)
#'
#'   - the data used in the EMMLi analysis (two elements - the subsampled data, and the corresponding models)
#'
#'   - the true level of subsampling
#'
#'   - (optional) the correlation matrix of the single best model.
#'
#' @export

subsampleEMMLi <- function(landmarks, fractions, models, min_landmark, aic_cut = 0,
                           return_corr = FALSE, nrep = NULL) {

  if (!is.null(nrep)) {
    if (length(fractions) > 1) {
      stop("Only one subsampling fraction is allowed when multiple subsampling simulations are run.")
    }
    fractions <- rep(fractions, nrep)
  }

  fitEmmli <- function(landmarks, models, fraction, min_landmark, aic_cut) {

    dat <- subsampleLandmarks(landmarks = landmarks, fraction = fraction, models = models,
                              min_landmark = min_landmark)
    c <- paleomorph::dotcorr(dat$landmarks)

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
