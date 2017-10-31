################################################################################
#' IRSAL
#'
#' @param atlas object of class "atlas" created by \link[Morpho]{createAtlas}
#' @param landmarks k x 3 x n array containing reference landmarks of the sample
#' or a matrix in case of only one target specimen.
#' @param initial_fixed The fixed points for the initial
#' \link[Morpho]{placePatch} analysis.
#' @param n The number of patched points to use in each iteration of IRSAL.
#' @param reps The number of times to repeat the IRSAL procedure.
#' @param write.out If TRUE then the patched points from each iteration, for
#' each species, will be written to the working directory. This will NOT include
#' the temporary anchor points.
#' @param sp If write.out is true, then this determines which species is
#' written out. Defaults to 1 (the first species in the dataset). 2 would be
#' the second species in the dataset etc.
#' @param ... Additonal arguments passed to \link[Morpho]{placePatch}
#' @name IRSAL
#' @export

IRSAL <- function(atlas, landmarks, initial_fixed, n, reps, write.out = FALSE,
                  sp = NULL, ...) {
  # Get the initial starting point.
  start <- original <- Morpho::placePatch(atlas = atlas,
                                          dat.array = landmarks,
                                          keep.fix = initial_fixed,
                                          ...)

  pb <- txtProgressBar(min = 0, max = reps, initial = 0)
  # Make a series of percentage increases to go through.
  pcs <- round(seq.int(from = 0.10, to = 0.90, length.out = reps), 2)

  res <- array(0, dim = dim(start))
  # Larger loop to loop over specimens.
  for (j in seq_len(dim(landmarks)[3])) {
    print(paste0("Specimen ", j))
    for (i in seq(reps)) {

      if (i == 1) {
        # Extract the patch points to sample from, from the initial patching.
        patch_idx <- (dim(landmarks)[1] + 1):dim(original)[1]
        pps <- start[patch_idx, , j]
        # calculate initial bending energy.

      } else {
        # Otherwise get the patch points to sample from from the patch from the
        # previous iteration.
        patch_idx <- (nrow(new_landmarks) + 1):nrow(t)
        pps <- t[patch_idx, ]
      }

      samples <- sort(sample(1:dim(pps)[1], round(dim(pps)[1] * pcs[i], 0)))

      # Take the relevant parts out of the atlas...
      a_patch <- atlas[["patch"]]
      a_lms <- atlas[["landmarks"]]
      a_fixed <- atlas[["keep.fix"]]

      # Set up the new info. Pathces stays the same.
      new_patch <- a_patch
      # Landmarks have the sampled patches added to them.
      # Landmarks should now be their original length + n.
      new_lms <- rbind(a_lms, a_patch[samples, ])
      # And the fixed points are defined
      new_fixed <- c(a_fixed, (nrow(a_lms) + 1):(nrow(a_lms) + length(samples)))

      # Now set up the new atlas. It's different because it has some new
      # fixed points (patches moved to landmarks)
      new_atlas <- atlas
      new_atlas[["patch"]] <- new_patch
      new_atlas[["landmarks"]] <- new_lms
      new_atlas[["keep.fix"]] <- new_fixed

      # Now add the new 10 datapoints onto the original data...
      new_landmarks <- abind::abind(landmarks[,,j], pps[samples, ], along = 1)
      prefix <- dimnames(landmarks)[3][[1]][j]
      # Now place the new patches.
      # function version
      # t <- Morpho::placePatch(atlas = new_atlas,
      #                         dat.array = new_landmarks,
      #                         keep.fix = new_fixed,
      #                         prefix = prefix,
      #                         ...)

      # testing version
      t <- Morpho::placePatch(atlas = new_atlas,
                              dat.array = new_landmarks,
                              keep.fix = new_fixed,
                              path="./ply",
                              prefix = prefix,
                              fileext = ".ply",
                              ray = TRUE,
                              inflate=2,
                              tol=5,
                              relax.patch = FALSE,
                              rhotol = NULL,
                              silent = FALSE,
                              mc.cores = 1)

      # calculate bending energy matrix...
      # This is from inside slider3d - find out what stepsize and m are...
      # This needs to be compared to the bending energy of the previous iteration
      # in some way...
      # L <- CreateL(t, output = "Lsubk3")
      # dataslido <- calcGamma(U$Gamma0, L$Lsubk3, U$U,
      #                        dims = m, stepsize = stepsize)

      if (write.out) {
        f_lms <- (nrow(landmarks[,,j]) + 1):(nrow(landmarks[,,j]) +
                                               length(samples))
        tt <- t[-f_lms, ]
        filename <- paste0("iteration_", i, "_", prefix, ".csv")
        write.csv(tt, file = filename, row.names = FALSE)
      }
      setTxtProgressBar(pb, i)
    }
    # Remove anchor points (patched points that are now landmarks) and store
    # result.
    f_lms <- (nrow(landmarks[,,j]) + 1):(nrow(landmarks[,,j]) + length(samples))
    t <- t[-f_lms, ]
    res[, , j] <- t
  }
  return(res)
}
