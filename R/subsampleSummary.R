###########################################################################################
#' getBestMods
#'
#' Collates the best fitting models from each of the analyses in a subsampledEMMLi analysis
#' @name getBestMods
#' @param rs The output of subsampleEMMLi
#' @keywords internal

getBestMods <- function(rs) {
  bestMods <- vector(mode = "list", length = length(rs))
  for (i in seq_along(rs)) {
    bestMods[[i]] <- cbind(t(rs[[i]][[1]]),
                           subsample = as.numeric(names(rs)[i]),
                           true_subsample = as.numeric(rs[[i]]$true_subsample))
  }
  bestMods <- do.call(rbind, bestMods)
  bestMods <- cbind(bestMods, subsample = as.numeric(names(rs)))
  return(bestMods)
}

###########################################################################################
#' getRhos
#'
#' Collates the rhos for each of the analyses in a subsampleEMMLi output. Collects them
#' together based on the model that derived them, like with like.
#' @name getRos
#' @param rs The output of subsampleEMMLi
#' @keywords internal

getRhos <- function(rs) {
  rho_names <- unlist(lapply(rs, function(x) names(x[[2]])))
  allRhos <- vector(mode = "list", length = length(unique(rho_names)))
  names(allRhos) <- unique(rho_names)
  for (i in seq_along(unique(rho_names))) {
    allRhos[[i]] <- vector(mode = "list", length = sum(rho_names == names(allRhos)[[i]]))
  }

  ref_names <- vector(mode = "list", length = length(allRhos))
  for (i in seq_along(rs)) {
    c_rs <- rs[[i]]$rhos_best
    for (j in seq_along(c_rs)) {
      slt <- which(names(allRhos) == names(c_rs)[j])
      first_null <- which(sapply(allRhos[[slt]], is.null))[[1]]
      xx <- cbind(c_rs[[j]],
                  subsample = as.numeric(names(rs)[i]),
                  true_subsample = as.numeric(rs[[i]]$true_subsample))
      bets <- grep("to", colnames(xx))
      ordered <- sapply(bets, function(x)
        paste(sort(as.numeric(strsplit(colnames(xx)[x], " to ")[[1]])), collapse = " to ")
      )
      # Sort out the right ordering for column names here...
      colnames(xx)[bets] <- ordered
      xx <- xx[ , order(colnames(xx), decreasing = TRUE)]
      xx <- xx[ , c(4:ncol(xx), 1:3)]
      allRhos[[slt]][[first_null]] <- xx
    }
  }

  for (i in seq_along(allRhos)) {
    allRhos[[i]] <- do.call(rbind, allRhos[[i]])
    allRhos[[i]] <- allRhos[[i]][rownames(allRhos[[i]]) != "MaxL", ]
  }
  return(allRhos)
}

###########################################################################################
#' subsampleSummary
#'
#' Summarises the output of subsampleEMMLi to allow easier comparison between multiple
#' subsampled analyses
#' @name summariseResults
#' @param subsamples The output of subsampleEMMLi
#' @return
#' @export

summariseResults <- function(subsamples) {
  bestMods <- getBestMods(subsamples)
  bestRhos <- getRhos(subsamples)
  return(list(bestModels = bestMods, bestRho = bestRhos))
}
