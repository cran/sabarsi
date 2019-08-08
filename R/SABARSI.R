#' Perform background removal on the whole SERS spectrum data set.
#' Divide the SERS spectrum data into time-frequency blocks and remove background locally.
#'
#' @param x A p*n data matrix. There are n SERS spectra with p frequency channels.
#' @param q A number taking value between 0 and 1. 100*q is the quantile that SABARSI uses to calculate the spectrum strength.
#' @param w.chan The window size for the frequency domain. The default value of \code{w.chan} is 50.
#' @param w.time The window size for the time domain. The default value of \code{w.time} is 50.
#' @return A p*n data matrix, \code{xr}, of background-removed spectra.
#' @import stats
#' @examples
#' data(SERS)
#' x <- SERS$R1   ## x is the matrix of SERS spectra
#' xr <- background_removal(x) ## xr is the matrix of background removed spectra
#' @export
background_removal <- function(x, q=0.4, w.chan=50, w.time=50)
{
  cent <- get.centers(nrow(x), ncol(x), w.chan, w.time)
  message(nrow(cent$st.grid), "blocks to calculate...")

  xr <- matrix(max(x) * 100, nrow=nrow(x), ncol=ncol(x))

  for (k in 1 : nrow(cent$st.grid))
  {
    st.chan <- cent$st.grid[k, 1]
    st.time <- cent$st.grid[k, 2]

    rgx <- (st.chan) : (st.chan + w.chan * 2 - 1)
    rgy <- (st.time) : (st.time + w.time * 2 - 1)

    ## get the residual of this calculation
    r.block <- norm.a.block(x[rgx, rgy], q)

    ## if the residual is larger, do not update. this keeps the results conservative.
    ## Note that the index operations here are quite tricky
    which.keep <- (abs(xr[rgx, rgy]) < abs(r.block))
    r.block[which.keep] <- (xr[rgx, rgy])[which.keep]
    xr[rgx, rgy] <- r.block
  }

  return(xr)
}

#' Calculate the pvalues and false discovery rates (FDRs) for a background-removed spectrum
#'
#' @param xr A p-length vector that represents a background-removed spectrum with p frequency channels.
#' @import stats
#' @return A list containing a vector of pvalues, \code{pvals}, and a vector of FDRs, \code{fdrs}.
#' @export

detect_sig <- function(xr)
{
  mad <- median(abs(xr)) / qnorm(0.75)
  pvals <- 2 * pnorm(-abs(xr) / mad)
  pvals[xr < 0] <- 1
  fdrs <- p.adjust(pvals, method = 'BH')
  return(list(pvals = pvals, fdrs = fdrs))
}

#' Detect signals in background-removed spectra
#'
#' @param xr A p*n data matrix of background-removed spectra, where n is the number of spectra, and p is the number of frequency channels.
#' @param fdr.c A number between 0 and 1, which is the cutoff for FDR.
#' @param stren.c A positive number for the cutoff of the strength of a signal.
#' @param wid.c A positive number for the cutoff of the width of a signal.
#' @import stats
#' @return  A list containing the indices of spectra with detected signals, \code{tim.index}, and a matrix recording the peaks of signals, \code{peaks}.
#' @examples
#' \donttest{
#' res <- signal_detection(xr)
#' # xr is the matrix of background removed spectra
#' # xr can be obtained by background_removal function
#'
#' head(res$tim.index)  ## check the first few time indices of signals
#' }
#' @export
signal_detection <- function(xr, fdr.c = 0.01, stren.c = 200, wid.c = 5)
{

  peaks <- data.frame(time=NA, st=NA, ed=NA, len=NA, minfdr=NA, pk=NA, cen=NA)
  n <- ncol(xr)
  p <- nrow(xr)
  fdr <- matrix(NA, nrow = p, ncol = n)

  ## calculate the fdrs for each background-removed spectrum
  for (i in 1:n) {
    fdr[,i] <- detect_sig(xr[,i])$fdrs
  }
  res <- group.sig.tol(xr,fdr,tol=fdr.c, wid.c = wid.c)
  peaks <- res$peaks
  res2 <- group.sig.tol(xr,fdr,tol=fdr.c, stren.cutoff = stren.c, wid.c = wid.c)
  tim.index <- res2$tim.index

  return(list(tim.index = tim.index, peaks = peaks, fdr = fdr, fdr.c = fdr.c))
}


#' Obtain the set of signature signals
#' Merge groups of concatenated signals and give the time indices of signature signals.
#'
#' @param xr A p-length vector that represents a background-removed spectrum with p frequency channels.
#' @param object An object obtained from \code{signal_detection}
#' @param t.tol A positive integer, which is the tolerance of time difference when comparing two signals.
#' @param cor.tol A number between 0 and 1. Two signals is considered to be similar enough if their correlation is higher than \code{cor.tol}.
#' @import stats
#' @return A vector recording the time indices of signature signals.
#' @examples
#' \donttest{
#' ## xr is the matrix of background removed spectra.
#' res <- signal_detection(xr) ## detect the signals in xr
#' tim.index.ss <- merge_signals(xr = xr, object = res)  ## the set of signature signals.
#' }
#' @export
#'
merge_signals<- function(xr, object, t.tol = 4, cor.tol = 0.5)
{
  tim.index <- object$tim.index
  n <- length(tim.index)
  peaks <- object$peaks
  fdr <- object$fdr
  fdr.c <- object$fdr.c
  xs <- simplify.spectrum(xr,peaks, tim.index)
  xss <- as.matrix(xr[,tim.index[1]], nrow=nrow(xr))
  tim.index.ss <- tim.index[1]
  for(i in 2:n)
  {
    if(tim.index[i] > (tim.index[i-1]+t.tol))
    {
      tim.index.ss <- c(tim.index.ss,tim.index[i])
      xss <- cbind(xss,xr[,tim.index[i]])
    }else if(corr_clean(xra = xr[,tim.index[i-1]], xrb = xr[, tim.index[i]], xs[,i-1], xs[,i]) > cor.tol)
    {
      sa <- median(xr[which(fdr[,tim.index[i-1]] < fdr.c), tim.index[i-1]])
      sb <- median(xr[which(fdr[,tim.index[i]] < fdr.c), tim.index[i]])
      if(sa < sb)
      {
        k<- length(tim.index.ss)
        tim.index.ss[k] <- tim.index[i]
        xss[,k] <- xr[,tim.index[i]]
      }
    }
    else{
      tim.index.ss <- c(tim.index.ss, tim.index[i])
      xss <- cbind(xss, xr[,tim.index[i]])
    }
  }
  return(tim.index.ss)
}

#' Match signals from two experiments.
#' For each signal in the first experiment, \code{shift.match} function finds the best matched signal in the second experiment. This function takes the potential frequency shifts into consideration for similarity measurement.
#'
#' @param xra A p*n data matrix of background-removed spectra in the first experiment.
#' @param xrb A p*n data matrix of background-removed spectra in the second experiment.
#' @param ta A vector of time indices of signals in the first experiment.
#' @param tb A vector of time indices of signals in the second experiment.
#' @import stats
#' @return A list containing the time indices of signals in the first experiment, \code{ta}, and the time indices of corresponding time indices in the second experiment, as well as the correlation of each match pairs, \code{corra}.
#' @export

shift_match<- function(xra, xrb, ta, tb)
{
  tb.m <- rep(0, length = length(ta))
  corra <- rep(0,length = length(ta))
  p <- nrow(xra)
  for(i in 1:length(ta))
  {
    corrb <- shift <- rep(0, length = length(tb))
    for(j in 1:length(tb))
    {
      corr.tmp <- rep(0,length=21)
      for(k in -10:10)
      {
        corr.tmp[k+11] <- cor(xra[15:(p-15), ta[i]], xrb[(15+k):(p-15+k), tb[j]])
      }
      corrb[j] <- max(corr.tmp)
      shift[j] <- which.max(corr.tmp) - 11
    }
    corra[i] <- max(corrb)
    tb.m[i] <- tb[which.max(corrb)]
  }
  return(data.frame(ta = ta, tb = tb.m, corelation = corra))
}

