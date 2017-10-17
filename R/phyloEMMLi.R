###########################################################################################
#' phyloEmmli
#' Takes landmarks and a phylogeny and then corrects the landmarks for the phylogeny according
#' to one of two methods, and then either returns the corrected landmarks, or calculates the
#' correlation matrix (using dotcorr) and fits EMMLi. Species missing from data or tree are
#' automatically dropped.
#' @name phyloEmmli
#' @param landmarks The landmarks. These can be in a 2D format with species as rownames,
#' and x, y, z as columns, or a 3D array with species names in the 3rd ([,,x]) dimension.
#' @param phylo A phylogeny describing the relationship between species in landmarks.
#' @param method Either "pgls" or "ic". If PGLS then corrected data are calculated as the
#' residuals of a phylogenetic least squares regression against 1 (\link[caper]{pgls}),
#' if IC then independent contrasts (\link[ape]{pic}).
#' @param EMMLi Logical - if TRUE EMMLi is fit and the results returned as well as the
#' phylogenetically correceted landmarks
#' @param ... Extra arguments required for EMMLi (at minimum models, and N_sample)
#' @export
#' @return A 2D or 3D array (depending on input) containing phylogenetically corrected
#' landmarks. If EMMLi = TRUE then the results of the EMMLi model are also returned.

phyloEmmli <- function(landmarks, phylo, method = "pgls", EMMLi = FALSE, ...) {
  # first check that the three is a tree.

  if (class(phylo) != "phylo") {
    stop("Tree must be an object of class 'phylo'.")
  }

  # now check that the data has rownames that are the species.
  if(length(dim(landmarks)) == 3) {
    if (is.null(dimnames(landmarks)[[3]])) {
      stop("Landmarks must have species names in the 3rd dimension.")
    } else {
      dims <- 3
    }
  } else if (length(dim(landmarks)) == 2) {
    if (is.null(rownames(landmarks))) {
      stop("Landmarks must have species names as rownames.")
    } else {
      dims <- 2
    }
  } else {
    stop("Landmarks must be either a 2D or 3D array.")
  }

  # check if mod and N_sample are provided if EMMLi is TRUE.
  if (EMMLi) {
    if (!exists(mod)) {
      stop("mod must be provided to fit EMMLi to phylo landmarks.")
    }
    if (!exists(N_sample)) {
      stop("N_sample must be provided to fit EMMLi to phylo landmarks.")
    }
  }

  # now check that the species on the tree match the species - here is where to check what format
  # the input is in, and check that species names are there.
  if (dims == 3) {
    sp_lm <- dimnames(landmarks)[[3]]
  } else if (dims == 2) {
    sp_lm <- rownames(landmarks)
  }

  sp_tr <- phylo$tip.label

  if (sum(sp_lm %in% sp_tr) != length(sp_lm)) {
    missing <- sum(!sp_lm %in% sp_tr)
    if (missing == length(sp_lm)) {
      stop("No species in dataset found on tree.")
    }
    print("Dropping", missing, "species from dataset - not in tree.")
    landmarks <- landmarks[sp_lm %in% sp_tr, ]
  }

  if (sum(sp_tr %in% sp_lm) != length(phylo$tip.label)) {
    missing <- sum(!sp_tr %in% sp_lm)
    if (missing == length(phylo$tip.label)) {
      stop("No species on tree found in dataset.")
    }
    print("Dropping", missing, "species from tree - not in dataset.")
    xtips <- phylo$tip.label[!sp_tr %in% sp_lm]
    phylo <- ape::drop.tip(phylo, xtips)
  }

  if (method == "pgls") {
    if (dims == 2) {
      landmarks <- as.data.frame(landmarks)
      landmarks$names <- rownames(landmarks)
      comp_data <- caper::comparative.data(phylo, landmarks, names = names)
      lms <- head(colnames(landmarks), -1)
      parallel::makeCluster(parallel::detectCores() - 2)
      parallel::clusterExport(cl, varlist = c("comp_data"), envir = environment())
      x <- parallel::parLapply(cl, lms, function(x) caper::pgls(formula(paste(x, "~", 1)), comp_data)$phyres)
      parallel::stopCluster(cl)
      phy_landmarks <- do.call(cbind, x)
      rownames(phy_landmarks) <- rownames(landmarks)
    } else if (dims == 3) {
      phy_landmarks <- landmarks
      for (i in seq_len(dim(landmarks)[2])) {
        m_lms <- as.data.frame(t(landmarks[,i,]))
        m_lms$names <- rownames(m_lms)
        comp_data <- caper::comparative.data(phylo, m_lms, names = names)
        lms <- head(colnames(m_lms), -1)
        cl <- parallel::makeCluster(parallel::detectCores() - 2)
        parallel::clusterExport(cl, varlist = c("comp_data"), envir = environment())
        x <- parallel::parLapply(cl, lms, function(x) caper::pgls(formula(paste(x, "~", 1)), comp_data)$phyres)
        parallel::stopCluster(cl)
        tmp_landmarks <- do.call(cbind, x)
        rownames(tmp_landmarks) <- rownames(m_lms)
        phy_landmarks[,i,] <- t(tmp_landmarks)
      }
    }
  } else if (method == "ic") {
    if (dims == 2) {
      landmarks <- landmarks[match(phylo$tip.label, rownames(landmarks)), ]
      phy_landmarks <- apply(landmarks, 2, function(x) ape::pic(x, phylo))
    } else (dims == 3) {
      phy_landmarks <- landmarks
      for (i in seq_len(dim(landmarks)[2])) {

      }
    }
  }

  # Fit or don't fit EMMLi.
  if (!EMMLi) {
    res <- phy_landmarks
  } else if (EMMLi) {
    # Here, arrayspecs may not be needed if the data is provided in 3D format.
    if (dims == 2) {
      arr <- geomorph::arrayspecs(phy_landmarks, ncol(phy_landmarks) / 3, 3)
    }
    corr <- paleomorph::dotcorr(arr)
    N_sample <- dim(landmarks)[3]
    emm <- EMMLi(...)
    res <- list(EMMLi = emm, phy_landmarks = phy_landmarks)
  }

  return(res)
}
