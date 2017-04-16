
library(FactoMineR)
library(DT)

heb <- read.table(file ="data-raw/hebergement.txt" )
colnames(heb) <- c("hotel", "location", "camping", "residencesec")
rownames(heb) <- c("agriculteur", "cadresup", "inactif", "ouvrier")

t <- vector()
for(i in 1:4)
  t <- append(t, sum(heb[,i]))
heb[5,] <- t
rownames(heb)[5] <- "total"

l <- vector()
for(i in 1:5)
  l <- append(l, sum(heb[i,]))
heb[,5] <- l
colnames(heb)[5] <- "total"

h <- heb[c(1:4), c(1:4)]
res <- CA(h)

total <- function(data){
  x <- colSums(data)
  data <- rbind(data, x)
  y <- apply(data, 1, sum)
  data <- cbind(data, y)
  colnames(data)[ncol(data)] <- 'total'
  rownames(data)[nrow(data)] <- 'total'
  return(data)
}

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


gcol <- function(data){
  data <- total(data)
  x <- data[nrow(data),1:ncol(data)-1]/data[nrow(data),ncol(data)]
  return(x)
}


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


valeurp <- function(data){
  k <- min(nrow(data), ncol(data))
  v <- ACP1(data)$values[2:k]
  return(v)
}


vecteurp <- function(data){
  k <- min(nrow(data), ncol(data))
  v <- ACP1(data)$vectors[,2:k]
  return(v)
}


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


plotPL <- function(data){
  c <- coProLigne(data)
  d <- coProColonne(data)
  plot(c[,1], c[,2], xlim = c(1.1*min(c[,1]), 1.1*max(c[,1])), ylim =c(1.1*min(c[,2]), 1.1*max(c[,2])),pch=19, col = "yellow", xlab = "", ylab = "", cex = 3)
  abline(v=0, col="black", lty=5)
  abline(h=0, col="black", lty=5)
  text(c[,1], c[,2], row.names(c), pos = c(1,2), col = "#006633", cex = 2)
}


plotPC <- function(data){
  d <- coProColonne(data)
  plot(d[,1], d[,2], xlim = c(1.1*min(d[,1]), 1.1*max(d[,1])), ylim =c(1.1*min(d[,2]), 1.1*max(d[,2])),pch=17, col = "red", xlab = "", ylab = "", cex = 3)
  abline(v=0, col="black", lty=5)
  abline(h=0, col="black", lty=5)
  text(d[,1], d[,2], row.names(d), cex = 2, pos = c(1,2), col = "#330066")
}



function(input, output, session) {
  output$plot <- renderPlot({
    if (input$plot_type == "Types of Accommodations") {
      plotPC(h)
      title(main = "Projection de variables sur plan de 2 dimensions", 
            xlab = "1er composant princial", ylab = "2eme composant princial",
            cex.main = 2,   font.main= 4
      )
    } else if (input$plot_type == "Socio-Professional Categories") {
      plotPL(h)
      title(main = "Projection des individus sur plan de 2 dimensions", 
            xlab = "1er composant princial", ylab = "2eme composant princial",
            cex.main = 2,   font.main= 4
      )
    } else if(input$plot_type == "Comparison") {
      plotAFC(h)
      title(main = "Le type de logement de vacance selon les CSP", 
            xlab = "1er composant princial", ylab = "2eme composant princial",
            cex.main = 2,   font.main= 4
      )
    }
    })

  output$summary <- renderPrint({
    summary(res)
  })

  output$table <- DT::renderDataTable({
    DT::datatable(heb)
  })
}
