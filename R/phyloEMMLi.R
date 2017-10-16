###########################################################################################
#' phyloEmmli
#' Takes landmarks and a phylogeny and then corrects the landmarks for the phylogeny according
#' to one of two methods, and then either returns the corrected landmarks, or calculates the
#' correlation matrix (using dotcorr) and fits EMMLi. Species missing from data or tree are
#' automatically dropped.
#' @name phyloEmmli
#' @param landmarks The landmarks. These can be in a 2D format with species as rownames, and x, y, z as columns, or a 3D array with species names in the 3rd ([,,x]) dimension.
#' @param phylo A phylogeny describing the relationship between species in landmarks.
#' @param method Either "pgls" or "ic". If PGLS then corrected data are calculated as the residuals of a phylogenetic least squares regression against 1, if IC then independent contrasts.
#' @param EMMLi Logical - if TRUE EMMLi is fit and the results returned.
#' @param ... Extra arguments required for EMMLi (at minimum models, and N_sample)
#' @export

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
  } else if (dim(landmarks) == 2) {
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
    # This bit needs changing for the case where landmarks are 3D.

    # When this is 3D I think I need to rebuild the data to a 2D form - where the
    # rows are the species names, and the columns are (xyz) per landmark, and then
    # work the phylogenetic correction over those columns. Otherwise for each species
    # I have a matrix... OR... well isn't it the same to just go over each combinations
    # of dimensions for the third, and get all of the data and do it that way?
    # I think that actually, yes, this is the same... so I just need to change the
    # way that the function is applied.

    lms <- colnames(landmarks)
    landmarks$names <- rownames(landmarks)
    comp_data <- caper::comparative.data(phylo, as.data.frame(landmarks), names = names)
    x <- lapply(lms, function(x) caper::pgls(formula(paste(x, "~", 1)), comp_data)$phyres)
    phy_landmarks <- do.call(cbind, x)
    rownames(phy_landmarks) <- rownames(landmarks)
  } else if (method == "ic") {
    # match the order of the dataset to the order of the tree.
    landmarks <- landmarks[match(phylo$tip.label, rownames(landmarks)), ]
    phy_landmarks <- apply(landmarks, 2, function(x) ape::pic(x, phylo))
  }

  # Fit or don't fit EMMLi.
  if (!EMMLi) {
    res <- phy_landmarks
  } else if (EMMLi) {
    # Here, arrayspecs may not be needed if the data is provided in 3D format.
    arr <- geomorph::arrayspecs(phy_landmarks, ncol(phy_landmarks) / 3, 3)
    corr <- paleomorph::dotcorr(arr)
    N_sample <- dim(landmarks)[3]
    emm <- EMMLi(...)
    res <- list(EMMLi = emm, phy_landmarks = phy_landmarks)
  }

  return(res)
}
