#' @export
Tabcet <- function(df){
  gc <- colSums(df)/nrow(df)
  sous <- function(x, y){return (x-y)}
  xc <- apply(df,1, sous, y=gc)
  xc <- t(xc)
  return(xc)
}

#' @export
Tabvarc <- function(df){
  xc <- Tabcet(df)
  xc <- as.matrix(xc)
  v <- (1/nrow(df))*t(xc)%*%xc
  return(v)
}

#' @export
Tabcetr <- function(df){
  v <- Tabvarc(df)
  xc <- Tabcet(df)
  d1v <- diag(sqrt(1/diag(as.matrix(v))))
  xr <- xc %*%d1v
  return(xr)
}

#' @export
Tabcorr <- function(df){
  v <- Tabvarc(df)
  d1v <- diag(sqrt(1/diag(as.matrix(v))))
  r <- d1v %*% v %*% d1v
  return(r)
}

#' @export
VALpro <- function(df){
  v <- Tabvarc(df)
  r <- eigen(v)
  return(r$values)
}

#' @export
VECTprop <- function(df){
  v <- Tabvarc(df)
  r <- eigen(v)
  return(r$vectors)
}

#' @export
compoPcl<- function(df){
  r <- Tabcorr(df)
  xr <- Tabcetr(df)
  c <- xr %*% eigen(r)$vectors
  return(c)
}

#' @export
cercleVA <- function(df){
  R <- Tabcorr(df)
  Rc1<-round(sqrt(eigen(R)$values[1])*(eigen(R)$vectors[,1]),3)
  Rc2<-round(sqrt(eigen(R)$values[2])*(eigen(R)$vectors[,2]),3)

  cX1<-c(Rc1[1],Rc2[1]) # coordonnées de X1
  cX2<-c(Rc1[2],Rc2[2]) # coordonnées de X2

  par(pty="s")
  plot(Rc1,Rc2,xlim=c(-1,1),ylim=c(-1,1),xlab="C1",ylab="C2")
  symbols(0,0,circles=1, inches=F, add=T)
  abline(v=c(0,0), h=c(0,0))
}

#' @export
coordVA <- function(df){
  R <- Tabcorr(df)
  Rc1<-round(sqrt(eigen(R)$values[1])*(eigen(R)$vectors[,1]),3)
  for (i in 2:length(eigen(R)$values)){
  Rc2<-round(sqrt(eigen(R)$values[i])*(eigen(R)$vectors[,i]),3)
  Rc1 <- cbind(Rc1, Rc2)}
  return(Rc1)
}

#' @export
cosinusVA <- function(df){
  c <- coordVA(df)
  d <- c^2
  return(d)
}

#' @export
contribuVA <- function(df){
  VECTprop(df) -> b
  l <- VALpro(df)
  l0 <- sqrt(l)
  prod <- function(x,y){return(x*y)}
  apply(b, 1, prod, l0) -> corr
  corr2 <- corr^2
  div <- function(x,y){return(x/y)}
  apply(100*corr2, 2, div, l) -> corr
  return(corr)
}

#' @export
coordIND <- function(df){
  compoPcl(df) -> l
  return(l)
}

#' @export
cosinusIND <- function(df){
  c <- coordIND(df)
  s <- rowSums(c^2)
  c2 <- c^2
  div <- function(x,y){return(x/y)}
  c <- apply(c2, 2, div, s)
  return(c)
}

#' @export
contribuIND <- function(df){
  c <- coordIND(df)
  l <- VALpro(df)
  c2 <- (1/nrow(df)) * (c^2)
  div <- function(x,y){return(x/y)}
  c <- apply(c2, 1, div, l)
  c <- t(c)
  return(c)
}


