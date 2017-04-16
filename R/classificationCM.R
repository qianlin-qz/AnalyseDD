#'@export
classCM <- function(data, k, g0){
  if (missing(g0)== TRUE) {
  g0 <- data[sample(1:nrow(data),k),]
  }
  g0 <- cbind(g0, "class" = 1:2)
  g0 <- g0[order(g0[,3]),]
  d0 <- distance3(data, g0)
  p0 <- pos(d0)
  c0 <- ltdf(p0)
  g1 <- ncentreg(data, p0, k)
  d1 <- distance3(data, g1)
  p1 <- pos(d1)
  c1 <- ltdf(p1)
  while (any(c1!=c0)){
    p0 <- p1
    c0 <- ltdf(p0)
    g2 <- ncentreg(data, p0, k)
    d2 <- distance3(data, g2)
    p1 <- pos(d2)
    c1 <- ltdf(p1);
    if (all(c1) == all(c0)) break;
    return(p1)
  }
  p1 <- ltdf(p1)
  p1 <- cbind(data, p1)
  p1 <- p1[order(p1[,3]),]
  colnames(p1) <- c(colnames(data), "class")
  return(list("classification" = p1, 'point de d¨¦part' = g0));
}

#'@export
# nouveau centre de gravit¨¦
ncentreg <- function(data, positn, k){
   r1 <- which(positn == 1)
   x1 <- data[r1, 1]
   y1 <- data[r1, 2]
   g <- c(sum(x1)/length(x1), sum(y1)/length(y1))
   for (i in k){
     r2 <- which(positn == i)
     gi <- c(sum(data[r2, 1])/length(data[r2, 1]), sum(data[r2, 2])/length(data[r2, 2]))
     g <- rbind(g, gi)
   }
   return(g)
}

#data frame to one element list
dftl <- function(data){
l <- list()
d <- for (i in 1:nrow(data)){
for (j in 1:ncol(data)){
l <- c(l,data[i,j])}}
return(l)
}

#'@export
# position de data en fonction de distance
pos <- function(data){
  data <- as.matrix(data)
  pst <- list()
  for (i in 1:nrow(data)){
    pst <- c(pst,which.min(data[i,]))
  }
  pst <- as.matrix(pst, ncol = 1)
  rownames(pst) <- 1:nrow(data)
  return(pst)
}

#'@export
# une table de tous les distances entre la data origine et le centre gravit¨¦
distance3 <- function(df1, df2){
  g1 <- c(df2[1,1], df2[1,2])
  c <- distance2(df1, g1)
  for (i in 2:nrow(df2)){
    gk <- c(df2[i,1], df2[i,2])
    d <- distance2(df1, gk)
    c <- cbind(c, d)
  }
  return(c)
}

#' @export
# distance entre deux point avec cordonn¨¦e (x,y)
distance2 <- function(df1, g){
  l1 <- lapply(g[[1]], distance1, y = df1[,1])
  l2 <- lapply(g[[2]], distance1, y = df1[,2])
  l <- cbind(as.data.frame(l1), as.data.frame((l2)))
  l <- as.data.frame(rowSums(l))
  return(l)
}

#'@export
#distance entre x et y en 1 dimension
distance1 <- function(x, y){
  return((x-y)^2)
}

#' @export
# list to dataframe
ltdf <- function(lst){
  lst <- as.data.frame(as.matrix(lst))
  lst <- as.list(lst)
  lst <- lapply(lst, unlist)
  lst <- as.data.frame(lst)
  return(lst)
}

