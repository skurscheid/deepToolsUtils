#fileParsers.R

# makeHTSeqCountMatrix()
# takes a vectory of HTSeq count output file names (including path)
# and returns a matrix
makeHTSeqCountMatrix <- function(files = NULL){
  if(is.null(files)) {stop("Filenames missing")}
  n <- length(files)
  for (i in 1:n){
    s <- unlist(lapply(strsplit(names(files)[i], "\\."), function(x) x[1]))
    if (i == 1) {
      mat0 <- read.table(files[i], header = F, as.is = T, sep = ("\t"))
      rn <- mat0[,i]
      mat0 <- as.matrix(mat0[,-1])
      colnames(mat0) <- s
      rownames(mat0) <- rn
    } else {
      mat1 <- read.table(files[i], header = F, as.is = T, sep = ("\t"))
      mat1 <- as.matrix(mat1[,-1])
      colnames(mat1) <- s
      mat0 <- cbind(mat0, mat1)
    }
  }
  return(mat0)
}
