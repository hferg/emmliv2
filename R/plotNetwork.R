###########################################################################################
#' plotNetwork
#'
#' Plots a network diagram of the output of an EMMLi analysis. Nodes are proportionally
#' sized to within-module rho, and lines are larger and darker according to between-module
#' rho.
#' @name plotNetwork
#' @param rhos The rhos that come out of an EMMLi analysis
#' @param module_names The names of the modules - if absent generic numbers are used.
#' If "rhos" then the nodes are named with the within-module rho.
#' @param linecolour The colour of the joining lines. If a single colour, then lines will
#' be that colour with width and transparency adjusted to the strength of the correlation.
#' If "viridis" then lines are coloured according to the viridis colour palette, with darker
#' colours corresponding to stronger correlations.
#' @param title Title for the plot.
#' @param layout A matrix describing the positions of each module on the canvas.
#' See qgraph::qgraph for details.
#' @return A plotted network of the relationships between and within modules.
#' @examples
#' emm <- EMMLi(corr = corr, mods = mods, N_sample = 34)
#' plotNetwork(emm$rho[[1]], linecolour = "viridis")
#' @export

plotNetwork <- function(rhos, module_names = NULL, linecolour = "#56B4E9",
  title = NULL, layout = NULL) {

  withins <- grep("Module*", colnames(rhos))
  nmodule <- length((grep("Module*", colnames(rhos))))
  rholist <- t(rhos)

  words <- strsplit(rownames(rholist), " ")

  plotcorr <- matrix(data = NA, nrow = nmodule, ncol = nmodule)
  modnums <- unlist(lapply(withins, function(x) strsplit(colnames(rhos)[x], " ")[[1]][2]))

  mods <- sapply(words[withins], function(x) paste(x, collapse = " "))
  mods <- gsub("Module", "M", mods)
  mods <- mods[order(sapply(mods, function(x)
    as.numeric(strsplit(x, "M ")[[1]][[2]])))]
  colnames(plotcorr) <- rownames(plotcorr) <- mods

  for (i in 1:(length(words) - 1)) {
    if (length(words[[i]]) == 2) {
      module <- paste("M", words[[i]][2])
      plotcorr[module, module] <- rholist[i, "MaxL_p"]
    }

    if (length(words[[i]]) == 3) {
      from_module <- paste("M", (words[[i]][1]))
      to_module <- paste("M", words[[i]][3])
      plotcorr[from_module, to_module] <- rholist[i, "MaxL_p"]
      plotcorr[to_module, from_module] <- rholist[i, "MaxL_p"]
    }
  }

  within <- diag(plotcorr)
  between <- plotcorr

  if (is.null(module_names)) {
    mod.names <- mods
  } else if (module_names == "rhos") {
    mod.names <- within
  }

  if (linecolour == "viridis") {
    cls <- viridis::viridis(100, direction = -1)
    vcols <- apply(between * 100, 1, function(x) cls[x])
    linecolour <- NULL
  } else {
    vcols <- NULL
  }

  qgraph::qgraph(between,
         shape = "circle",
         posCol = linecolour,
         edge.color = vcols,
         labels = mod.names,
         vsize = within * 10,
         diag = FALSE,
         title = title,
         layout = layout)
}

###########################################################################################
#' plotRandomSubsamples
#'
#' Generates network plots for a given number of randomly selected subsampled EMMLi analyses,
#' as returned by subSampleEmmli.
#' @name plotRandomSubsamples
#' @param subsasmples The output of the execution of subSampleEmmli
#' @param n The number of random subsamples to plot.
#' @param ... Additional arguments for plotNetwork (see ?plotNetwork)
#' @return A gridded plot with a network plot for each of n random subsamples.
#' @examples
#' ssemm <- subsampleEmmli(landmarks = landmarks, models = models, fractions = 0.4, min_landmark = 5, nsim = 25)
#' plotRandomSubsamples(subsamples = ssemm, n = 9, linecolour = "viridis")
#' @export

plotRandomSubsamples <- function(subsamples, n, ...) {
  samples <- sample(1:length(subsamples$results), n)
  par(mfrow = n2mfrow(n))
  for (i in samples) {
    rh <- subsamples$results[[i]]$rhos_best[[1]]
    plotNetwork(rh, ...)
  }
}
