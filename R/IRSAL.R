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
#' @param write_out If TRUE then the landmarks for the first species each
#' iteration will be written to the working directory.
#' @param sp If write.out is true, then this determines which species is
#' written out. Defaults to 1 (the first species in the dataset). 2 would be
#' the second species in the dataset etc.
#' @param ... Additonal arguments passed to \link[Morpho]{placePatch}
#' @export

#' New info: When estimated patch points are moved into landmarks they should
#' stay there (I think?) and in such a way more points move over to fixed as
#' the analysis progresses. Andre says he starts with 15% and ends up at around
#' 90% - which implies that over 200 iterations some attempts must be rejected.
#' It looks like he uses an estimate of bending energy to decide whether to
#' keep or reject a subsample - that can be calculated during the sliding
#' step of which normally takes place after this.
#' the function that looks like it will get me there is slider3d
#'

# Take some points from the patch, add them to landmarks.
# Project new patch using the new landmarks.
# If bending energy is improved, take a slightly larger number of points, and
# add them to the original landmarks (i.e. NOT the landmarks from the previous
# step).
# Repeat, increasing the number sampled slightly each time. It is this
# increasing step that is important.


IRSAL <- function(atlas, landmarks, initial_fixed, n, reps, write_out = FALSE,
                  sp = NULL, ...) {
  # Get the initial starting point.
  start <- original <- Morpho::placePatch(atlas = atlas,
                                          dat.array = landmarks,
                                          keep.fix = initial_fixed,
                                          ...)
  lms_dim <- dim(landmarks)
  org_dim <- dim(original)
  patch_idx <- (lms_dim[1] + 1):org_dim[1]
  pb <- txtProgressBar(min = 0, max = reps, initial = 0)
  # Make a series of percentage increases to go through.
  pcs <- round(seq.int(from = 0.15, to = 0.90, length.out = reps), 2)
  for (i in seq(reps)) {

    # Get patch points to sample from.
    # If first iteration, get them from start.
    if (i == 1) {
      # Extract the patch points to sample from, from the initial patching.
      pps <- start[patch_idx, , ]
      # calculate initial bending energy.

    } else {
      # Otherwise get the patch points to sample from from the patch from the
      # previous iteration.
      # I THINK that the ordering is not necesary in light of chatting to Andre
      # since the patches are just added back onto the initial landmarks...
      # First order according to patch_ref.
      # reorder t.
      t <- t[match(ref_order$og, ref_order$new_order), , ]
      pps <- t[patch_idx, , ]
    }

    samples <- sort(sample(1:dim(pps)[1], n))

    # define the reference order for putting the "fixed" patches back where
    # they belong.
    ref_order <- data.frame(
      og = c(1:org_dim[[1]]),
      new_order = c(
        c(1:lms_dim[[1]]),
        (samples + lms_dim[[1]]),
        c(patch_idx[!patch_idx %in% (samples + lms_dim[[1]])]))
    )

    a_patch <- atlas[["patch"]]
    a_lms <- atlas[["landmarks"]]
    a_fixed <- atlas[["keep.fix"]]

    new_patch <- a_patch[samples * -1, ]
    new_lms <- rbind(a_lms, a_patch[samples, ])
    new_fixed <- c(a_fixed, (nrow(a_lms) + 1):(nrow(a_lms) + n))

    new_atlas <- atlas
    new_atlas[["patch"]] <- new_patch
    new_atlas[["landmarks"]] <- new_lms
    new_atlas[["keep.fix"]] <- new_fixed

    # Now add the new 10 datapoints onto the original data...
    new_landmarks <- abind::abind(landmarks, pps[samples, , ], along = 1)

    t <- Morpho::placePatch(atlas = new_atlas,
                            dat.array = new_landmarks,
                            keep.fix = new_fixed,
                            ...)

    if (write_out) {
      if (is.null(sp)) {
        sp <- 1
      } else {
        sp <- sp
      }
      t <- t[match(ref_order$og, ref_order$new_order), , ]
      write.csv(t[,,sp], file = paste0("iteraion", i, ".csv"),
                row.names = FALSE)
    }
    setTxtProgressBar(pb, i)
  }
  t <- t[match(ref_order$og, ref_order$new_order), , ]
  return(t)
}
