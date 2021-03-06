% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subsampleEMMLi.R
\name{subsampleEMMLi}
\alias{subsampleEMMLi}
\title{subsampleEMMLi}
\usage{
subsampleEMMLi(landmarks, fractions, models, min_landmark, aic_cut = 0,
  return_corr = FALSE, nrep = NULL)
}
\arguments{
\item{landmarks}{A 2D array of xyz landmarks to subsample from. Will be turned into a correlation
matrix using \link[paleopmorph]{dotcorr} for EMMLi analysis after subsampling.}

\item{fractions}{Either a single subsampling fraction (in which case nrep is required) or a
vector of fractions. Specified in decimal format, i.e. a fraction of 0.2 will subsample down
to 20\% of the original number of landmarks.}

\item{models}{A data frame defining the models. The first column should contain the landmark names
as factor or character. Subsequent columns should define which landmarks are contained within each
module with integers, factors or characters. If a landmark should be ignored for a specific model
(i.e., it is unintegrated in any module), the element should be NA.}

\item{min_landmark}{The minimum number of landmarks to subsample to. This ensures that a module isn't
totally removed during random subsampling. When subsampling causes a landmark to be removed or to be
subsampled below this threshold landmarks are drawn from the original module at random and added back in.
This means that sometimes (especially with low subsampling fractions and/or the presence of small
modules in a model) the actual subsampling level is higher than the requested subsampling. In these
cases a warning is printed to the screen.}

\item{aic_cut}{This is the threshold of dAICc below which two models are considered to be not different.
When this occurs multiple models are be returned as the best. Defaults to 0 (i.e., only the best
fitting model is returned, regardless of how close other models are.)}

\item{return_corr}{Logical - if TRUE then the full correlations of within and between modules are
returned for the single best model in addition to the EMMLi results. Defaults to FALSE.}

\item{nrep}{If a single subsampling fraction, this is the number of times that subsampling fraction
is repeated.}
}
\value{
A list of n elements, where n is either the number of subsampling fractions, or nrep. Each
element contains the output of an EMMLi analysis on the subasampled data, consisting of four or five
elements:

  - the best model(s) description

  - the rho list(s) for the best model(s)

  - the data used in the EMMLi analysis (two elements - the subsampled data, and the corresponding models)

  - the true level of subsampling

  - (optional) the correlation matrix of the single best model.
}
\description{
Analyse random subsamples of a dataset repeatedly using EMMLi. This function subsamples
a 2D array of landmarks to a given fraction of it's original size and then fits EMMLi
to the smaller dataset, returning the results. This can be either done repeatedly for
a single subsampling fraction, or for a range of subsampling fractions (e.g. to investigate
the effects of increasing subsampling). In all cases EMMLi is fit to a correlation matrix
calculated from the subsampled landmarks using \link[paleomorph]{dotcorr}.
}
\details{
If a single fraction is provided then an nrep argument is also required, and the function
will subsample the data to the given fraction nrep times. Alternatively, if a range of
fractions is given (in a vector) the nrep argument is not required, and the function will
subsample the data and fit EMMLi once for each fraction in the given fractions vector.
}
