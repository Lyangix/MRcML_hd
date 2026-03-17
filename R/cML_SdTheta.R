#' Standard Error of Estimated Theta
#'
#' Get the standard error of estimated theta using RAPS-style variance
#' (mr.raps.simple formula) for better calibration. Uses only valid IVs when K>0.
#'
#' @param b_exp Vector of estimated effects for exposure.
#' @param b_out Vector of estimated effects for outcome.
#' @param se_exp Vector of standard errors for exposure.
#' @param se_out Vector of standard errors for outcome.
#' @param theta Estimated theta from cML.
#' @param b_vec Estimated vector of b from cML (unused; kept for compatibility).
#' @param r_vec Estimated vector of r from cML (identifies valid IVs when r_vec==0).
#'
#' @return Standard error of theta.
#' @keywords internal
cML_SdTheta <- function(b_exp, b_out,
                        se_exp, se_out,
                        theta, b_vec, r_vec)
{
  zero_ind <- which(r_vec == 0)
  if (length(zero_ind) < 2L) return(NaN)

  b_exp <- b_exp[zero_ind]
  b_out <- b_out[zero_ind]
  se_exp <- se_exp[zero_ind]
  se_out <- se_out[zero_ind]

  denom <- (se_out^2 + theta^2 * se_exp^2)^2
  score_var <- sum(((b_exp^2 - se_exp^2) * se_out^2 + (b_out^2 - se_out^2) * se_exp^2 + se_exp^2 * se_out^2) / denom)
  I <- sum(((b_exp^2 - se_exp^2) * se_out^2 + (b_out^2 - se_out^2) * se_exp^2) / denom)

  if (I <= 0 || score_var < 0) return(NaN)
  sqrt(score_var / I^2)
}
