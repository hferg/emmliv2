###########################################################################################
#' compareModules
#'
#' This takes a correlation matrix or 3D landmark array, a model definition, and then two module
#' numbers or names to compare. It plots (if plot = TRUE) a figure of three boxplots - the first
#' two are the correltions within each of the two modules, and the third is the between-module
#' correlations. These boxes are coloured such that matching colours are not significantly
#' different according to a Tukey HSD test. The results of the anova and tukey HSD test are
#' also returned.
#' @name compareModules
#' @param corr A correlation matrix or a 3D array of landmarks. If 3D then a correlation matrix
#' is calculated with \link[paleomorph]{dotcorr}
#' @param model Either a vector of numbers describing a model of modules, or a 2 column dataframe
#' with the first bein landmark names and the second being the module definitions.
#' @param test_modules A vector of two module numbers to compare, or if the modules are named, the names
#' of those two modules.
#' @param plot Logical - if TRUE the plot is drawn.
#' @return A list with two elements - the first is the result of an ANOVA compaing the mean
#' correlations within- and between-modules, and the second is the results of a TukeyHSD test on
#' that ANOVA. If plot = TRUE a plot showing these results is called.
#' @export
#' @examples
#' data(macacaCorrel)
#' data(macacaModels)
#' # Pick a model to draw modules from - as a vector.
#' model <- macacaModels$Goswami
#' compareModules(corr = macacaCorrel, model = model, test_modules = c(2, 5))
#'
#' # Or as a 2 column dataframe...
#' model <- macacaModels[ , c(1, 4)]
#' compareModules(corr = macacaCorrel, model = model, test_modules = c(2, 5))

compareModules <- function(corr, model, test_modules, plot = TRUE) {

  # If landmarks supplied
  if (length(dim(corr)) == 3) {
    require(paleomorph)
    corr <- paleomorph::dotcorr(corr)
  }

  # If model has landmark names
  if (is.vector(model)) {
    lms <- model
  } else if (ncol(model) == 2) {
    lms <- array(model[ , 2])
    # If model is just a vector
  } else if (ncol(model) > 2) {
    stop("Model must either be a vector of model definitions or a data frame with the first column containing landmark names and the second of module definitions.")
  }

  symmet = corr
  symmet[upper.tri(symmet)] = t(symmet)[upper.tri(symmet)]

  if (is.vector(model)) {
    modNF <- stats::na.omit(cbind(1:length(model), lms))
  } else {
    modNF <- stats::na.omit(cbind(1:nrow(model), lms))
  }

  w <- unique(modNF[, 2])
  w <- w[w %in% test_modules]
  all_modules <- list()
  modules <- list()
  btw_mod = list()
  betweenModules = list()
  withinModules = list()
  unintegrated = list()
  betweenFloat = list()
  for(i in seq(length(w))){
    fg <- modNF[modNF[, 2] == w[i], ]
    l <- corr[as.numeric(fg[, 1]), as.numeric(fg[, 1])]
    modules[[i]] <- (as.array(l[!is.na(l)]))
  }

  names(modules) <- paste("Module", w)
  if (length(w) > 1) {
    cb <- utils::combn(w, 2)
    for (i in seq(dim(cb)[2])){
      fg1 <- modNF[modNF[, 2] == cb[1, i], ]
      fg2 <- modNF[modNF[, 2] == cb[2, i], ]
      btw <- symmet[
        as.integer(setdiff(fg2[, 1], fg1[, 1])),
        as.integer(setdiff(fg1[, 1], fg2[, 1]))]
      btw_mod[[i]] <- btw[!is.na(btw)]
    }
    names(btw_mod) <- paste(cb[1, ], "to", cb[2, ])
    betweenModules['betweenModules'] = list(as.vector(rle(unlist(btw_mod))$values))
    withinModules['withinModules'] = list(as.vector(rle(unlist(modules))$values))
  }
  all_modules = c(modules, btw_mod, withinModules, betweenModules)

  # Prepare data.
  td <- c(all_modules[[1]], all_modules[[2]], all_modules[[3]])
  group <- as.factor(c(
    rep(names(all_modules)[1], length(all_modules[[1]])),
    rep(names(all_modules)[2], length(all_modules[[2]])),
    rep(names(all_modules)[3], length(all_modules[[3]]))
  ))
  td <- data.frame(corrs = td, group = group)
  a <- aov(corrs ~ group, data = td)
  t <- TukeyHSD(a)

  res <- list(anova = a, tukeyhsd = t)

  if (plot) {
    # function to group variables that are not different.
    t_labs <- function(t, v){
      levs <- t[[v]][,4]
      labs <- data.frame(multcompView::multcompLetters(levs)['Letters'])
      labs$treatment <- rownames(labs)
      labs <- labs[order(labs$treatment) , ]
      return(labs)
    }
    labels <- t_labs(t, "group")

    # add labels to td for colouring.
    td$color <- NA
    for (i in seq_len(nrow(labels))) {
      td$color[td$group == labels$treatment[i]] <- as.character(labels$Letters[i])
    }

    pall <- c("#E69F00", "#56B4E9", "#7BB31A")
    p <- ggplot2::ggplot(td, aes(x = group, y = corrs, fill = color)) +
      geom_boxplot() +
      scale_fill_manual(values = pall) +
      xlab("") +
      ylab("Correlation") +
      guides(fill = guide_legend(title = "Significance")) +
      theme(
        legend.position = "none"
      )+
      stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                   width = .75, linetype = "dashed")
    print(p)
  }
  return(res)
}
