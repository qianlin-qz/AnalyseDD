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

write.table(heb, file = "data/hebergement.csv", append = FALSE, row.names = TRUE)
save(heb,  file = "data/hebergement.rdata")
