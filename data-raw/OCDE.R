read.csv(file = "data-raw/OCDE.txt", sep = "", dec = ",") -> ocde

rownames(ocde) <- ocde[,2]
div <- function(x,y){return(x/y)}
ocde <- apply(ocde[, 4:12], 2, div, ocde[,3])
ocde <- round(ocde, 4)


save(ocde,  file = "data/OCDE.rdata")
