#' MRcMLhd: Mendelian Randomization via Constrained Maximum Likelihood (High-Dimensional)
#'
#' Implements the high-dimensional version of the MRcML method for Mendelian
#' Randomization. The key difference from the standard \pkg{MRcML} package is
#' the use of a RAPS-style sandwich variance estimator in \code{cML_SdTheta},
#' which provides better-calibrated standard errors when the number of genetic
#' instruments is large.
#'
#' The two main user-facing functions are:
#' \itemize{
#'   \item \code{\link{mr_cML_hd}} — cML without data perturbation.
#'   \item \code{\link{mr_cML_DP_hd}} — cML with data perturbation.
#' }
#'
#' @keywords internal
"_PACKAGE"
