# functions to load, summarise and visualise computeMatrix output
require(jsonlite)
require(readr)
require(GenomicRanges)
require(IRanges)



computeMatrixLoader <- function(matrix_file = NULL){
  stopifnot(!is.null(matrix_file))
  if (file.exists(matrix_file)){
    # load first line of matrix, as it contains JSON run defininition
    runDef <- readLines(gzfile(matrix_file),n = 1)
    runDef <- jsonlite::fromJSON(gsub("@", "", runDef))

    computeMatrix <- readr::read_tsv(gzfile(matrix_file), comment = "@", col_names = FALSE)
    gr <- GenomicRanges::GRanges(seqnames = computeMatrix$X1,
                                 IRanges::IRanges(start = computeMatrix$X2,
                                                  end = computeMatrix$X3,
                                                  names = computeMatrix$X4),
                                 strand = computeMatrix$X6)
    computeMatrix <- computeMatrix[,-c(1:3,5:6)]
    fn <- unlist(strsplit(matrix_file, "/"))[length(unlist(strsplit(matrix_file, "/")))]
    runList <- list(runDef = runDef,
                    computeMatrix = computeMatrix,
                    gr = gr,
                    fileName = fn)
    return(runList)
  } else {
    stop("File does not exist")
  }
}

covPlotStats <- function(x){
  if (is.list(x)){
    x <- unlist(x)
  }
  l <- length(x)
  x <- as.numeric(x)
  m <- mean(x)
  stdev <- sd(x)
  sem <- sd(x)/sqrt(l)
  ci_lower <- m-2*sem
  ci_upper <- m+2*sem
  rl <- data.frame(length = l,
             mean = m,
             stdev = stdev,
             sem = sem,
             ci_lower = ci_lower,
             ci_upper = ci_upper)
  return(rl)
}

makePlottingData <- function(computeMatrixList = NULL){
  stopifnot(!is.null(computeMatrixList))
  stopifnot(is.list(computeMatrixList))
  cm <- computeMatrixList
  bins <- (cm$runDef$upstream + cm$runDef$downstream) / cm$runDef$`bin size`
  maxBin <- max(cm$runDef$sample_boundaries)
  l1 <- lapply(cm$runDef$sample_boundaries, function(x){
    start <- x + 2
    end <- start + bins - 1
    if (x < maxBin){
      st <- cm$computeMatrix[,c(start:end)]
      return(st)
    }
  })
  l1 <- l1[1:length(cm$runDef$sample_labels)]
  names(l1) <- cm$runDef$sample_labels
  l2 <- lapply(names(l1), function(x){
    print(x)
    sumDat <- apply(l1[[x]], 2, covPlotStats)
    sumDat <- do.call("rbind", sumDat)
    m1 <- data.frame(bin = c(1:bins),
                     value = sumDat$mean,
                     stdev = sumDat$stdev,
                     sem = sumDat$sem,
                     ci_lower = sumDat$ci_lower,
                     ci_upper = sumDat$ci_upper,
                     n_genes = sumDat$length,
                     sample = x,
                     group = paste(unlist(strsplit(x, "-"))[1:2], collapse = "-"))
    return(m1)
  })
  m1 <- as.data.frame(do.call("rbind", l2))
  return(m1)
}



