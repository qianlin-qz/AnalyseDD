#'Totaux
#'
#'Ajouter la colonne et la ligne totaux
#'@param data: data frame construite par importer les data ¨¤ traiter
#'@return un tableau de data d'origine avec les totaux lignes et colonnes
#'@export
total <- function(data){
  x <- colSums(data)
    data <- rbind(data, x)
    y <- apply(data, 1, sum)
    data <- cbind(data, y)
    colnames(data)[ncol(data)] <- 'total'
    rownames(data)[nrow(data)] <- 'total'
  return(data)
}

#'Tableau des profils-lignes
#'
#'Le tableau des profils en lignes est construit en divisant l'effectif de chaque case par le total de la ligne correspondante
#'@param data: valeur d'origine
#'@return un tableau des profils-lignes
#'@export
tabProLigne <- function(data){
    data <- total(data)
    x <- data.frame(matrix(NA, nrow = nrow(data), ncol=ncol(data)))
    for (i in 1:ncol(data)){
      x[,i] <- data[,i]/data[,ncol(data)]
    }
    colnames(x) <- colnames(data)
    rownames(x) <- rownames(data)
    x <- round(x, 3)
    return(x)
  }

#'Tableau des profils-colonnes
#'
#'Le tableau des profils en colonnes est construit en divisant l'effectif de chaque case par le total de la colonne correspondante
#'@param data: valeur d'origine
#'@return un tableau des profils-colonnes
#'@export
tabProColonne <- function(data){
    data <- total(data)
    x <- data.frame(matrix(NA, nrow = nrow(data), ncol=ncol(data)))
    for (i in 1:nrow(data)){
      x[i,] <- data[i,]/data[nrow(data),]
    }
    colnames(x) <- colnames(data)
    rownames(x) <- rownames(data)
    x <- round(x, 3)
    return(x)
  }

#'Tableau des valeurs theoriques
#'
#'Le tableau des valeurs th¨¦oriques est construit en divisant le total de ligne fois le total de colonne par le total de tableau
#'@param data: data frame construite par importer les data ¨¤ traiter
#'@return un tableau des valeurs th¨¦oriques
#'@export
tabValTheo <- function(data){
    data <- total(data)
    x <- data.frame(matrix(NA, nrow = nrow(data), ncol=ncol(data)))
    for (i in 1:ncol(data)){
      x[,i] <- data[nrow(data),i]*data[, ncol(data)]/data[nrow(data),ncol(data)]
    }
    colnames(x) <- colnames(data)
    rownames(x) <- rownames(data)
    return(x)
  }

#' Le test d'independance
#'
#' Pour mesurer de la liasion entre deux variables qualitatives, on realiser un test khi-2. Cette statistique est la r¨¦alisation d'une variable al¨¦atoire d2 de la loi khi-2 avec degre de libert¨¦ (r-1)(c-1)
#' @param data : valeur d'origine
#' @return la valeur de d2, la valeur critique et la conclusion de test
#' @export
d2 <- function(data){
  vo <- total(data)
  vt <- tabValTheo(data)
  val <- sum((vt-vo)^2/vt)
  valc <- qchisq(0.95, (nrow(data)-1)*(ncol(data)-1))
  if (val > valc){
    return(list("X-squared" = val,"valcritique" = valc, "Summary" = "Les variables sont dependantes"))
  }
  else {
    return(list("X-squared" = val,"valcritique" = valc, "Summary" = "Les variables sont independantes."))
  }
}

#'Tableau de contributions au khi-2
#'
#'@param data : valeur d'origine
#'@return contribution au khi-2 de chaque observation
#'@export
contrX2 <- function(data){
  vo <- total(data)
  vt <- tabValTheo(data)
  d2 <- d2(data)[[1]]
  x <- ((vt-vo)^2/vt)/d2
  x[nrow(x),] <- sapply(x[1:nrow(x)-1,],sum)
  for (i in 1:nrow(x)){
  x[i,ncol(x)] <- sum(x[i,1:ncol(x)-1])
  }
  x <- round(100*x,2)
  return(x)
}

#'Tableau de frequences
#'
#'Le tableau de fr¨¦quences est construit en divisant le nombre observ¨¦ de la variable par le nombre total d'observation
#'@param data: valeur d'origine
#'@return un tableau de fr¨¦quences
#'@export
tabFreq <- function(data){
  data <- total(data)
  x <- data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data)))
  for (i in 1:ncol(data)){
    x[,i] <- data[,i]/data[nrow(data),ncol(data)]
  }
  colnames(x) <- colnames(data)
  rownames(x) <- rownames(data)
  x <- round(x, 3)
  return(x)
}

#'Centre gravite de colonnes
#'
#'Le centre de gravite de profils-lignes affect¨¦s avec des poids(frequences)
#'@param data: valeur d'origine
#'@return un vector de centre de gravite
#'@export
gcol <- function(data){
  data <- total(data)
  x <- data[nrow(data),1:ncol(data)-1]/data[nrow(data),ncol(data)]
  return(x)
}

#'Centre gravite de lignes
#'
#'Le centre de gravite de profils-colonnes affect¨¦s avec des poids(frequences)
#'@param data: valeur d'origine
#'@return un vector de centre de gravite
#'@export
grow <- function(data){
  data <- total(data)
  x <- data[1:nrow(data)-1,ncol(data)]/data[nrow(data),ncol(data)]
  return(x)
}

ACP1 <- function(data){
  n <- sum(rowSums(data))
  M1half <- diag(sqrt(n/colSums(data)))
  P1 <- tabProLigne(data)[1:nrow(data), 1:ncol(data)]
  P1 <- as.matrix(P1)
  M1half <- as.matrix(M1half)
  X <- P1 %*% M1half
  ACP1 <- eigen(t(X) %*% diag(rowSums(data)/n) %*% X)
  return(ACP1)
}

#'Valeur propre
#'
#'@param data : valeur d'origine
#'@return valeur propre
#'@export
valeurp <- function(data){
  k <- min(nrow(data), ncol(data))
  v <- ACP1(data)$values[2:k]
  return(v)
}

#'Vecteur propre
#'
#'@param data : valeur d'origine
#'@export
vecteurp <- function(data){
  k <- min(nrow(data), ncol(data))
  v <- ACP1(data)$vectors[,2:k]
  return(v)
}

#'Coordonnees des profils-lignes
#'
#'@param data : valeur d'origine
#'@export
coProLigne <- function(data){
  n <- sum(rowSums(data))
  M1half <- diag(sqrt(n/colSums(data)))
  P1 <- tabProLigne(data)[1:nrow(data), 1:ncol(data)]
  P1 <- as.matrix(P1)
  M1half <- as.matrix(M1half)
  X <- P1 %*% M1half
  Y <- vecteurp(data)
  C1 <- X %*% Y
  return(C1)
}

#'Coordonnees des profils-colonnes
#'
#'@param data : valeur d'origine
#'@export
coProColonne <- function(data){
  C1 <- coProLigne(data)
  C1 <- as.matrix(C1)
  P2 <- tabProColonne(data)
  P2 <- P2[1:nrow(data), 1:ncol(data)]
  P2 <- as.matrix(P2)
  P2 <- t(P2)
  Y <- valeurp(data)
  C2 <- P2 %*% C1 %*% diag(1/sqrt(Y))
  return(C2)
}

#'Graphique en utilisant la methode AFC
#'
#'@param data : valeur d'origine
#'@export
plotAFC <- function(data){
  c <- coProLigne(data)
  d <- coProColonne(data)
  if (max(c[,1])>max(d[,1]) && max(c[,2])>max(d[,2])){
  plot(c[,1], c[,2], xlim = c(1.1*min(c[,1]), 1.1*max(c[,1])), ylim =c(1.1*min(c[,2]), 1.1*max(c[,2])),pch=19, col = "yellow",  xlab = "", ylab = "", cex = 3)
  abline(v=0, col="black", lty=5)
  abline(h=0, col="black", lty=5)
  points(d[,1], d[,2],pch = 17, col = "red", cex = 3)
  text(c[,1], c[,2], row.names(c), pos = c(1,2), col = "#006633", cex = 2)
  text(d[,1], d[,2], row.names(d), pos = c(1,2), col = "#330066", cex = 2)
  }
  else if (max(c[,1])<max(d[,1]) && max(c[,2])<max(d[,2])){
    plot(d[,1], d[,2], xlim = c(1.1*min(d[,1]), 1.1*max(d[,1])), ylim =c(1.1*min(d[,2]), 1.1*max(d[,2])),pch=19, col = "yellow",  xlab = "", ylab = "", cex = 3)
    abline(v=0, col="black", lty=5)
    abline(h=0, col="black", lty=5)
    points(c[,1], c[,2], pch = 17, col = "red")
    text(c[,1], c[,2], row.names(c), pos = c(1,2), col = "#006633", cex = 2)
    text(d[,1], d[,2], row.names(d), pos = c(1,2), col = "#330066", cex = 2)
  }
  else {
    par(mfrow = c(1,2))
    plot(c[,1], c[,2], xlim = c(1.1*min(c[,1]), 1.1*max(c[,1])), ylim =c(1.1*min(c[,2]), 1.1*max(c[,2])), pch =19, col = "yellow",  xlab = "", ylab = "", cex = 3)
    plot(d[,1], d[,2], xlim = c(1.1*min(d[,1]), 1.1*max(d[,1])), ylim =c(1.1*min(d[,2]), 1.1*max(d[,2])),pch = 17, col = "red",  xlab = "", ylab = "", cex = 3)
    text(c[,1], c[,2], row.names(c), pos = c(1,2), col = "#006633", cex = 2)
    text(d[,1], d[,2], row.names(d), pos = c(1,2), col = "#330066", cex = 2)
  }
}

#'@export
plotAAA <- function(data){
  c <- coProLigne(data)
  d <- coProColonne(data)
  plot(c[,1], c[,2], xlim = c(1.1*min(c[,1]), 1.1*max(c[,1])), ylim =c(1.1*min(c[,2]), 1.1*max(c[,2])),pch=1, col = "red", cex.lab = 0)
  text(c[,1], c[,2], row.names(c), cex = 0.5, pos = c(1,2), col = "#006633")
  abline(v=0, col="black", lty=5)
  abline(h=0, col="black", lty=5)
  par(new=TRUE)
  plot(d[,1], d[,2], xlim = c(1.1*min(d[,1]), 1.1*max(d[,1])), ylim =c(1.1*min(d[,2]), 1.1*max(d[,2])),pch=1, col = "red", xlab = "", ylab = "")
  text(d[,1], d[,2], row.names(d), cex = 0.5, pos = c(1,2), col = "#330066")
}

#'@export
plotPL <- function(data){
  c <- coProLigne(data)
  d <- coProColonne(data)
  plot(c[,1], c[,2], xlim = c(1.1*min(c[,1]), 1.1*max(c[,1])), ylim =c(1.1*min(c[,2]), 1.1*max(c[,2])),pch=19, col = "yellow", xlab = "", ylab = "", cex = 3)
  abline(v=0, col="black", lty=5)
  abline(h=0, col="black", lty=5)
  text(c[,1], c[,2], row.names(c), pos = c(1,2), col = "#006633", cex = 2)
}

#'@export
plotPC <- function(data){
  d <- coProColonne(data)
  plot(d[,1], d[,2], xlim = c(1.1*min(d[,1]), 1.1*max(d[,1])), ylim =c(1.1*min(d[,2]), 1.1*max(d[,2])),pch=17, col = "red", xlab = "", ylab = "", cex = 3)
  abline(v=0, col="black", lty=5)
  abline(h=0, col="black", lty=5)
  text(d[,1], d[,2], row.names(d), cex = 2, pos = c(1,2), col = "#330066")
}


