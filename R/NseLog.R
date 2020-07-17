# LogNSE
NseLog <- function (m, o, nullModel=mean(o, na.rm=na.rm), na.rm=TRUE) {
    m <- log(m + 1e-04)
    o <- log(o + 1e-04)
    err1 <- sum((m - o)^2, na.rm=na.rm)
    err2 <- sum((o - nullModel)^2, na.rm=na.rm)
    ns <- 1 - (err1/err2)
    ns
}

