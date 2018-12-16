ganova <-
function(y,X,A,xi) {
  n <- length(y); r <- qr(X)$rank
  proj <- function(M) return(M%*%ginv(M)) # orthog. projection matrix
  R02 <- sum((y - proj(X)%*%y)^2)
  errordf <- n - r
  Xhyp <- X - X%*%ginv(A)%*%A
  yhyp <- y - X%*%ginv(A)%*%xi
  RH2 <- sum((yhyp - proj(Xhyp)%*%yhyp)^2)
  totaldf <- n - qr(Xhyp)$rank
  hypdf <- totaldf - errordf
  RegSS <- RH2 - R02
  RMS <- R02 / errordf; RegMS <- RegSS / hypdf
  Fstat <- RegMS / RMS
  SS <- c(RegSS,R02,RH2)
  DF <- c(hypdf,errordf,totaldf)
  MS <- c(RegMS,RMS,NA)
  Fratio <- c(Fstat,NA,NA)
  pvalue <- c(pf(Fstat,hypdf,errordf,lower.tail=F),NA,NA)
  Source <- c("Departure from H0","Error","Total")
  return(data.frame(cbind(SS,DF,MS,Fratio,pvalue),
                    row.names = Source))
}
