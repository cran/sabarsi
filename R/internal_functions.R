## get centers
get.centers <- function(len.chan, len.time, w.chan, w.time)
{
  ## number of slots on the channels
  n.chan <- floor(len.chan / w.chan)
  if (len.chan == n.chan * w.chan)
  {
    n.chan <- n.chan - 1
  }

  ## number of slots on the time
  n.time <- floor(len.time / w.time)
  if (len.time == n.time * w.time)
  {
    n.time <- n.time - 1
  }

  ## construct the slots on the channels and times
  st.chan <- c((0 : (n.chan-2)) * w.chan + 1, len.chan - w.chan * 2 + 1)
  st.time <- c((0 : (n.time-2)) * w.time + 1, len.time - w.time * 2 + 1)

  st.grid <- expand.grid(st.chan, st.time)

  return(list(st.chan=st.chan, st.time=st.time, st.grid=st.grid))
}

## q is the quantile we want to use
norm.a.block <- function(x.b, q)
{
  ## first, normalize each time point to get the graph roughly on the same scale
  time.ave <- colMeans(x.b)
  x.b.time <- t(t(x.b) / time.ave)

  ## then, normalize at channels
  bg <- apply(x.b.time, 1, median)

  ## then, get the projection of original data on the bg
  x.b.proj <- (x.b / bg)

  ## use the q quantile to get the strength of background
  s.proj <- apply(x.b.proj, 2, quantile, q)

  ## get what's left
  sigl <- x.b - bg %*% t(s.proj)

  return(sigl)
}

get.st.ed <- function(fdrs, xrr, fdr.cutoff)
{
  st <- which((c(fdrs, 1) < fdr.cutoff) & (c(1, fdrs) >= fdr.cutoff))     ###   the start point of a signal
  ed <- which((c(1, fdrs) < fdr.cutoff) & (c(fdrs, 1) >= fdr.cutoff)) - 1  ### the end point of a signal

  if (length(st) != length(ed))
  {
    stop("Lengths of st and ed are not the same!")  ### infact the starting length should be the same as the end length
  }

  len <- ed - st + 1   ### a vector of the length of each bump
  pk <- minfdr <- cen <- rep(1, length(st))     #### length(st)=length(ed)=# of bumps




  if (length(ed) > 0)
  {
    for (i in 1 : length(st))
    {
      if (st[i] > ed[i])  ### st[i] <= ed[i]
      {
        stop("Some of st is larger than ed!")
      }
      minfdr[i] <- min(fdrs[st[i] : ed[i]])    ### the min fdr of each bump: how significant a bump can be
      pk[i] <- max(xrr[st[i] : ed[i]])         ### the strongest intensity of each bump: how strong a bump's signal can be
      cen[i] <- which.max(xrr[st[i] : ed[i]]) + st[i] - 1      ### the center pixel of each bump
    }
  }

  return(data.frame(st=st, ed=ed, len=len, minfdr=minfdr, pk=pk, cen=cen))   ### return the imformation of all bumps' starting pixel,

}

########## this function is to return the features with significant bumps >= 5

group.sig.tol <- function(xr,fdr,tol=0.01,stren.cutoff=0, wid.c = 5)
{
  peaks <- data.frame(time=NA, st=NA, ed=NA, len=NA, minfdr=NA, pk=NA, cen=NA)
  icur <- 0
  fdr[which(xr<stren.cutoff)] <- 1
  for (j in 1 : ncol(fdr))     #### j is the index of significant time point's signal
  {
    tmp <- get.st.ed(fdr[, j], xr[, j], tol)   ### for each signal's fdr
    if (nrow(tmp) > 0)    ### nrow= # of signal bumps. if at this time point has any significant signal bumps
    {
      rag <- (icur + 1) : (icur + nrow(tmp))
      peaks[rag, 2 : ncol(peaks)] <- tmp
      peaks$time[rag] <- j   ### which time point the signal bump belongs to
      icur <- icur + nrow(tmp)
    }
  }
  peaks <- peaks[1 : icur,]
  tim.index<- unique(peaks$time[which(peaks$len>= wid.c)])
  return(list(peaks=peaks,tim.index=tim.index))
}


corr_clean<- function(xra,xrb,xsa,xsb)
{
  xu<- xsa+xsb
  #plot(xu,type="l")
  xu[which(xu>0)]<-1
  xra[which(xra<0)]<- 0
  xrb[which(xrb<0)]<- 0
  xra<- xra*xu
  xrb<- xrb*xu

  return(cor(xra,xrb))
}


simplify.spectrum<- function(xr, peaks, tim.index)
{
  n<- length(tim.index)
  xs<- matrix(0,ncol=n,nrow=nrow(xr))
  peaks<- peaks[peaks$len>=5,]
  for(i in 1:n)
  {
    ww<- which(peaks$time==tim.index[i])
    for(j in ww)
    {
      xs[(peaks$st[j]:peaks$ed[j]),i]<- 1
    }
  }
  return(xs)
}
