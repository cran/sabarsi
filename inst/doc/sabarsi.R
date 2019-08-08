## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(sabarsi)

## ------------------------------------------------------------------------
data(SERS)
x <- list()
x[[1]] <- SERS$R1
x[[2]] <- SERS$R2
dim(x[[1]])

## ------------------------------------------------------------------------
t1 <- 80
t2 <- 63
layout(matrix(1:2,nrow=2))
par(mar=c(0,0,0,0))
plot(x[[1]][,t1], type = "l", xlab = "Frequency", ylab = "Intensity")
plot(x[[2]][,t2], type = "l", xlab = "Frequency", ylab = "Intensity")

## ------------------------------------------------------------------------
xr <- list()
for (i in 1:2) {
  xr[[i]] <- background_removal(x[[i]])
}

## ------------------------------------------------------------------------
layout(matrix(1:2, nrow = 2))
par(mar=c(0,0,0,0))
plot(xr[[1]][,t1], type = "l", xlab = "Frequency", ylab = "Intensity", main = "R1")
plot(xr[[2]][,t2], type = "l", xlab = "Frequency", ylab = "Intensity", main = "R2")


## ------------------------------------------------------------------------
res <- list()
for (i in 1:2) {
  res[[i]] <- signal_detection(xr[[i]])
}

## ------------------------------------------------------------------------
head(res[[1]]$tim.index)

## ----message=FALSE-------------------------------------------------------
tim.index.ss <- list()
for (i in 1:2) {
  tim.index.ss[[i]] <- merge_signals(xr = xr[[i]], object = res[[i]])
}


## ------------------------------------------------------------------------
plot(xr[[1]][,t1], type = "l", col = "red", xlab = "Frequency", ylab = "Intensity")
lines(xr[[2]][,t2], col = "blue")

print(cor(xr[[1]][, t1], xr[[2]][, t2]))

## ------------------------------------------------------------------------

res.match <- shift_match(xr[[1]], xr[[2]], tim.index.ss[[1]],tim.index.ss[[2]])
print(res.match)

