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
#' @importFrom magrittr %>%
#' @name IRSAL
#' @export

IRSAL <- function(atlas, landmarks, initial_fixed, n, reps, write.out = FALSE,
                  sp = NULL, ...) {

  # Place patches to get a starting point for each specimen.
  print("Initial patching...")
  start <- original <- Morpho::placePatch(atlas = atlas,
                                          dat.array = landmarks,
                                          keep.fix = initial_fixed,
                                          silent = TRUE,
                                          ...)

  # Make a series of percentage increases to go through.
  pcs <- round(seq.int(from = 0.10, to = 0.90, length.out = reps), 2)

  res <- array(0, dim = dim(start))
  dimnames(res) <- dimnames(start)

  for (j in seq_len(dim(landmarks)[3])) {
    prefix <- dimnames(landmarks)[[3]][j]
    print(paste0("Specimen ", j, ":", prefix))
    t <- start[,,j]
    for (i in seq(reps)) {

      if (i == 1) {
        # If it's the first iteration, calculate bending energy relative to
        # the template. Otherwise the bending energy from the previous iteration
        # is used.
        template <- rbind(atlas[["landmarks"]], atlas[["patch"]])
        L_int <- CreateL(template)
        be_p <- t(start[,,j]) %*% L_int$Lsubk %*% start[,,j] %>%
                as.matrix %>%
                diag %>%
                sum
      }

      # Since when bending energy is tested the extra landmarks are removed,
      # the original dimentions of the landmarks are sufficient to determine
      # where the patches are in a new t (i.e. they are always the same)
      patch_idx <- (dim(landmarks)[1] + 1):dim(original)[1]
      pps <- t[patch_idx, ]

      # Sample some of the patches.
      samples <- sort(sample(1:dim(pps)[1], round(dim(pps)[1] * pcs[i], 0)))

      # Take the relevant parts out of the original atlas...
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
      # fixed points (landmarks is longer - patches moved there).
      new_atlas <- atlas
      new_atlas[["patch"]] <- new_patch
      new_atlas[["landmarks"]] <- new_lms
      new_atlas[["keep.fix"]] <- new_fixed

      # Now move the patches into the landmarks data to make new landmark
      # data to correspond to the definition in the new atlas.
      new_landmarks <- abind::abind(landmarks[,,j], pps[samples, ], along = 1)
      # Now place the new patch
      t_new <- Morpho::placePatch(atlas = new_atlas,
                              dat.array = new_landmarks,
                              keep.fix = new_fixed,
                              prefix = prefix,
                              silent = TRUE,
                              ...)

      # Compare bending energy to previous iteration.
      # Remove the "IRSALed" landmarks...
      f_lms <- (nrow(landmarks[,,j]) + 1):(nrow(landmarks[,,j]) + length(samples))
      t_test <- t_new[-f_lms, ]

      # Then compare the bending energy. Then compare this to the bending
      # energy of the last iteration (be_p)
      L <- CreateL(t_test)
      be_n <- t(t) %*% L$Lsubk %*% t %>%
              as.matrix %>%
              diag %>%
              sum

      # Now, if be_n (bending energy new) is smaller than be_p (bending
      # energy previous), t becomes t_test (the new patches, with the false
      # anchor points removed), and be_p becomes be_n.
      if (be_n < be_p) {
        t <- t_test
        be_p <- be_n
      }

      if (write.out) {
        f_lms <- (nrow(landmarks[,,j]) + 1):(nrow(landmarks[,,j]) +
                                               length(samples))
        tt <- t[-f_lms, ]
        filename <- paste0("iteration_", i, "_", prefix, ".csv")
        write.csv(tt, file = filename, row.names = FALSE)
      }
    }
    # Remove anchor points (patched points that are now landmarks) and store
    # result. At this point t will be the patch that has the lowest bending
    # energy during IRSAL - and it will also have had the landmarks removed
    # alread (it will either be the initial patch that hasn't ever changed,
    # if bending energy was never reduced, OR it will be a t_test (that has
    # had the IRSAL landmarks removed in order to compare bending energy) that
    # had a lower bending energy. Pop that patch into the overall results...
    res[, , j] <- t
  }
  return(res)
}
