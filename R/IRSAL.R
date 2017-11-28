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
#' @param maxit The maximum number of attempts to find improved bending energy.
#' @param ... Additonal arguments passed to \link[Morpho]{placePatch}
#' @importFrom magrittr %>%
#' @name IRSAL
#' @export

IRSAL <- function(atlas, landmarks, initial_fixed, n, reps, write.out = FALSE,
                  sp = NULL, maxit = 50, ...) {

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

  # Calculate the bending energy for the template to reference against.
  template <- rbind(atlas[["landmarks"]], atlas[["patch"]])
  L_int <- CreateL(template)

  for (j in seq_len(dim(landmarks)[3])) {
    prefix <- dimnames(landmarks)[[3]][j]
    print(paste0("Specimen ", j, ":", prefix))
    t <- start[,,j]

    # calculate initial bending energy
    be_p <- t(t) %*% L_int$Lsubk %*% t %>%
    as.matrix %>%
      diag %>%
      sum
    failure <- 0
    i <- 1
    repeat {

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
      be_n <- t(t_test) %*% L_int$Lsubk %*% t_test %>%
              as.matrix %>%
              diag %>%
              sum

      # Now, if be_n (bending energy new) is smaller than be_p (bending
      # energy previous), t becomes t_test (the new patches, with the false
      # anchor points removed), and be_p becomes be_n.
      if (be_n < be_p) {
        print(paste0("Bending energy improved at rep ", i))
        t <- t_test
        be_p <- be_n
        # Increment i
        i <- i + 1
        # Then check if it equals reps + 1 (which means it has fully gone
        # through all the sample percentages).
        if (i == (reps + 1)) {
          break
        }

      } else {
        failure <- failure + 1

        if (failure == maxit) {
          print(paste0("No resolution at rep ", i, " after ", maxit,
                       " attempts."))
          break
        }

      }
    }
    # Remove anchor points (patched points that are now landmarks) and store
    # result. At this point t will be the patch that has the lowest bending
    # energy during IRSAL - and it will also have had the landmarks removed
    # alread (it will either be the initial patch that hasn't ever changed,
    # if bending energy was never reduced, OR it will be a t_test (that has
    # had the IRSAL landmarks removed in order to compare bending energy) that
    # had a lower bending energy. Pop that patch into the overall results...
    print("Storing best patch.")
    res[, , j] <- t
  }
  return(res)
}
