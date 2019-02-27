append_kt_model <- function(obj, n_foreward=30, spe=c("ar1", "ma1")) {
  kt <- obj$model$kt
  start <- as.numeric(names(kt)[1])
  kt <- ts(kt)
  kt <- c(kt, sim_kt(kt, n_foreward=n_foreward, spe))
  names(kt) <- start + (1:length(kt)) - 1
  kt <- ts(kt, start=start)

  obj$forecast <- list(
    ax=obj$model$ax,
    bx=obj$model$bx,
    kt=kt,
    kt.spe=spe
  )
  obj
}



identify_kt <- function(obj) {
  kt <- obj$model$kt
  times <- 1:length(kt)
  res_kt <- resid(lm(kt~times))
  par(mfrow=c(1,2))
  acf(res_kt)
  pacf(res_kt)
  par(mfrow=c(1,1))
  eacf(ts(kt), 5, 5)
}


sim_kt <- function(kt, n_foreward=10, spe=c("ar1", "i1", "ma1")) {
  spe <- match.arg(spe)
  kt <- ts(kt)

  if (spe == "ar1") {
    m <- forecast::Arima(kt, c(1, 0, 0), include.drift=T)
    co <- m$coef
    r0 <- kt[length(kt)] - length(kt) * co["drift"] - co["intercept"]

    sim <- stats::arima.sim(list(ar=c(co["ar1"])), n_foreward, sd=sqrt(m$sigma2), start.innov=r0, n.start = 1)
    sim <- sim + (length(kt) + 1:n_foreward) * co["drift"] + co["intercept"]

  } else if (spe=="ma1") {
    m <- forecast::Arima(kt, c(0, 0, 1), include.drift=T)
    co <- m$coef
    r0 <- m$residuals[length(kt)]

    sim <- stats::arima.sim(list(ma=c(co["ma1"])), n_foreward, sd=sqrt(m$sigma2), start.innov=r0, n.start = 1)
    sim <- sim + (length(kt) + 1:n_foreward) * co["drift"] + co["intercept"]
  } else {
    m <- forecast::Arima(kt, c(0, 1, 0), include.drift=T)

    sim <- stats::arima.sim(list(order=c(0, 1, 0)), n_foreward, sd=sqrt(m$sigma2))[-1]
    sim <- sim + (length(kt) + 1:n_foreward) * m$coef["drift"] + m$coef["intercept"]
  }

  sim
}


sim_boot_kt <- function(kt, n_foreward=10, B=100, spe=c("ar1", "i1", "ma1")) {
  spe <- match.arg(spe)
  xs <<- ts(kt)

  ys <- matrix(0, B, n_foreward)

  if (spe == "ar1") {
    m <- forecast::Arima(xs, c(1, 0, 0), include.drift=T)
    drift <- m$coef["drift"]
    r0 <- kt[length(kt)] - length(kt) * drift - m$coef["intercept"]

    boo <- TSA::arima.boot(m, B=B)
    for (i in 1:B) {
      co <- boo[i,]
      sim <- stats::arima.sim(list(ar=c(co["ar1"])), n_foreward, sd=sqrt(co["noise var"]), start.innov=r0, n.start = 1)
      sim <- sim + (length(kt) + 1:n_foreward) * drift + co["intercept"]
      ys[i, ] <- sim
    }
  } else if (spe=="ma1") {
    m <- forecast::Arima(kt, c(5, 0, 0), include.drift=T)
    co <- m$coef
    r0 <- tail(kt, 5) - (length(kt)-4:0) * drift - m$coef["intercept"]

    boo <- TSA::arima.boot(m, B=B)
    for (i in 1:B) {
      co <- boo[i,]
      sim <- stats::arima.sim(list(ar=c(co[1:5])), n_foreward, sd=sqrt(co["noise var"]), start.innov=r0, n.start = 5)
      sim <- sim + (length(kt) + 1:n_foreward) * drift + co["intercept"]
      ys[i, ] <- sim
    }
  } else {
    m <- forecast::Arima(kt, c(0, 1, 0), include.drift=T)

    for (i in 1:B) {
      sim <- stats::arima.sim(list(order=c(0, 1, 0)), n_foreward, sd=sqrt(m$sigma2))[-1]
      sim <- sim + (length(kt) + 1:n_foreward) * m$coef["drift"] + m$coef["intercept"]
      ys[i, ] <- sim
    }
  }

  colnames(ys) <- 1:n_foreward
  rownames(ys) <- 1:B
  ys
}
