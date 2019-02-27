#' Title
#'
#' @param obj
#' @param type
#'
#' @return
#' @export
#'
#' @examples
residuals.LCM <- function(obj, type=c("deviance", "pearson")) {
  type <- match.arg(type)
  dhat <- obj$model$fitted
  ddat <- obj$data$N/obj$data$P

  if (obj$link=="log") {
    if (type=="deviance") {
      rd <- sign(ddat-dhat)*sqrt(ddat*(log(ddat)-log(dhat)) - (ddat-dhat))
    } else {
      rd <- (log(ddat) - log(dhat))*sqrt(ddat)
    }
  } else {
    if (type=="deviance") {
      # todo
    } else {
      # todo
    }
  }

  rd[is.na(rd)] <- 0
  rd[is.infinite(rd)] <- 0
  dimnames(rd) <- dimnames(ddat)

  rd
}
