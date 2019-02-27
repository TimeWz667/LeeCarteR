

bootstrap <- function(obj, B=1000, n_foreward=10, B_kt=B/10, spe=c("ar1", "ma1", "i1"), calculate=T) {
  spe <- match.arg(spe)
  n_p <- ceiling(B / B_kt)

  num <- obj$data$Event
  pop <- obj$data$Population
  n_per <- nrow(pop)
  n_age <- ncol(pop)
  n_obj <- n_per * n_age


  bot <- lapply(1:n_p, function(x=0) {
    if (obj$link == "log") {
      samp <- matrix(rpois(n_obj, num), n_per, n_age)
      res <- fit_lcm_poisson(samp, pop, fitted=F, se=F, max_iter=200)
    } else {
      samp <- matrix(rbinom(n_obj, pop, num/pop), n_per, n_age)
      res <- fit_lcm_binomial(samp, pop, fitted=F, se=F, max_iter=200)
    }
    res$kt.fore <- sim_kt(res$kt, n_foreward, spe)
    res$kt.fore.mc <- sim_boot_kt(kt=res$kt, n_foreward=n_foreward, B=B_kt, spe)
    res
  })


  mc <- list()
  i = 1

  for (b in 1:n_p) {
    ax <- bot[[b]]$ax
    bx <- bot[[b]]$bx
    kt <- bot[[b]]$kt
    for(k in 1:B_kt) {
      kt.fore <- bot[[b]]$kt.fore.mc[k, ]
      names(kt.fore) <- tail(obj$data$Periods, 1) + 1:n_foreward

      mc[[i]] <- list(
        ax=ax, bx=bx, kt=kt, kt.fore=kt.fore
      )
      i <- i + 1
    }
  }

  if (calculate) {
    mc <- lapply(mc, function(xt){
      if (obj$link == "log") {
        mx <- with(xt, exp(t(ax+bx%*%t(kt))))
        mx.fore <- with(xt, exp(t(ax+bx%*%t(kt.fore))))
      } else {
        mx <- with(xt, {1/(1+exp(t(ax+bx%*%t(kt))))})
        mx.fore <- with(xt, 1/(1+exp(t(ax+bx%*%t(kt)))))
      }
      colnames(mx) <- colnames(mx.fore) <- obj$data$Ages
      c(xt, list(mx=mx, mx.fore=mx.fore))
    })
  }
  mc
}





