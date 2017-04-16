comp <- read.table(file ="data-raw/comparaison.txt" )

colnames(comp) <- c("Q1O", "Q1N", "Q2O", "Q2N")

save(comp,  file = "data/comparaison.rdata")
