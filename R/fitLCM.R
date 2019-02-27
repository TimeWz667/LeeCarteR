#' Title
#'
#' @param num
#' @param pop
#' @param link
#' @param se
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
fit_lcm <- function(num, pop, periods=NULL, ages=NULL, link=c("log", "logit"), se=T, ...) {
  link <- match.arg(link)

  pop[pop<num] <- num[pop<num] + 1

  n_age <- ncol(pop)
  n_per <- nrow(pop)

  if (is.null(ages)) {
    ages <- colnames(pop)
    if (is.null(ages)) ages <- 1:n_age
  }

  if (is.null(periods)) {
    periods <- suppressWarnings(as.numeric(rownames(pop)))
    if (any(is.na(periods))) periods <- 1:n_per
  }

  if (link == "log") {
    model <- fit_lcm_poisson(num, pop, se, ...)
  } else {
    model <- fit_lcm_binomial(num, pop, se, ...)
  }


  res <- list(
    model=model,
    link=link,
    data=list(Event=num, Population=pop, Ages=ages, Periods=periods)
  )

  class(res) <- "LCM"

  res
}


fit_lcm_poisson <- function(num, pop, se=T, fitted=T, tol=10e-5, max_iter=1000) {
  n_age <- ncol(pop)
  n_per <- nrow(pop)

  update_ax <- function(ax, bx, kt) {
    d_hat <- pop*exp(t(ax+bx%*%t(kt)))
    grad <- colSums(num-d_hat)/colSums(d_hat)
    grad[is.na(grad)] <- 0
    grad[is.infinite(grad)] <- 0
    ax + grad
  }


  update_kt <- function(ax, bx, kt) {
    d_hat <- pop*exp(t(ax+bx%*%t(kt)))
    tbx <- matrix(bx, n_per, n_age, byrow=T)
    grad <- rowSums((num-d_hat)*tbx) / rowSums(d_hat*tbx*tbx)
    grad[is.na(grad)] <- 0
    grad[is.infinite(grad)] <- 0
    grad[1] <- 0
    kt <- kt + grad
    kt
  }

  update_bx <- function(ax, bx, kt) {
    d_hat <- pop*exp(t(ax+bx%*%t(kt)))
    tkt <- matrix(kt, n_per, n_age)
    grad <- colSums((num-d_hat)*tkt)/colSums(d_hat*tkt*tkt)
    grad[is.na(grad)] <- 0
    grad[is.infinite(grad)] <- 0
    grad[1] <- 0
    bx <- bx + grad
    bx / sum(bx)
  }

  fn_log_li <- function(ax, bx, kt) {
    d_hat <- pop*exp(t(ax+bx%*%t(kt)))
    d_hat[d_hat==0] <- 1
    sum(dpois(num, d_hat, log=T), na.rm=T)
  }

  # Initialise
  ax0 <- colSums(num/pop, na.rm=T)
  bx0 <- rep(1, n_age) / n_age
  kt0 <- rep(0, n_per)
  li0 <- fn_log_li(ax0, bx0, kt0)
  i <- 0

  while(i < max_iter) {
    ax1 <- update_ax(ax0, bx0, kt0)
    kt1 <- update_kt(ax1, bx0, kt0)
    bx1 <- update_bx(ax1, bx0, kt1)

    li1 <- fn_log_li(ax1, bx1, kt1)
    ax0 <- ax1; bx0 <- bx1; kt0 <- kt1

    if (li1 <= li0*(1+(sign(li0))*tol)) break
    li0 <- li1
    i <- i + 1
  }

  # Normalise
  scl <- sum(bx1)
  bx <- bx1/scl
  kt <- kt1*scl

  sft <- mean(kt1)
  kt <- kt1 - sft
  ax <- ax1 + sft*bx

  names(ax) <- names(bx) <- colnames(pop)
  names(kt) <- rownames(pop)

  res <- list(ax=ax, bx=bx, kt=kt)

  if (se) {
    axli <- function(x) {
      fn_log_li(x, bx, kt)
    }

    bxli <- function(x) {
      x <- c(1-sum(x), x)
      fn_log_li(ax, x, kt)
    }

    ktli <- function(x) {
      x <- c(-sum(x), x)
      fn_log_li(ax, bx, x)
    }

    ax.se <- sqrt(abs(diag(solve(-hessian(axli, ax)))))
    bx.se <- sqrt(abs(diag(solve(-hessian(bxli, bx[-1])))))
    bx.se <- c(0, bx.se)
    kt.se <- sqrt(abs(diag(solve(-hessian(ktli, kt[-1])))))
    kt.se <- c(0, kt.se)

    names(ax.se) <- names(bx.se) <- colnames(pop)
    names(kt.se) <- rownames(pop)
    res <- c(res, list(ax.se=ax.se, bx.se=bx.se, kt.se=kt.se))
  }



  if (fitted) {
    res$fitted <- exp(t(ax+bx%*%t(kt)))
    dimnames(res$fitted) <- dimnames(pop)
  }
  res
}


fit_lcm_binomial <- function(num, pop, se=T, fitted=T, tol=10e-5, max_iter=1000) {

}
