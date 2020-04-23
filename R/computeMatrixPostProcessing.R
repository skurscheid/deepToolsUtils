# functions to load, summarise and visualise computeMatrix output
require(jsonlite)
require(readr)
require(GenomicRanges)
require(IRanges)
require(data.table)

computeMatrixLoader <- function(matrix_file = NULL){
  stopifnot(!is.null(matrix_file))
  if (file.exists(matrix_file)){
    # load first line of matrix, as it contains JSON run defininition
    runDef <- readLines(gzfile(matrix_file),n = 1)
    runDef <- jsonlite::fromJSON(gsub("@", "", runDef))

    computeMatrix <- data.table::fread(matrix_file, skip = "1", sep = "\t", header = FALSE)
    gr <- GenomicRanges::makeGRangesFromDataFrame(computeMatrix, keep.extra.columns = T, seqnames.field = "V1", start.field = "V2", end.field = "V3", strand.field = "V6")
    rows <- computeMatrix[,4]
    computeMatrix <- computeMatrix[,-c(1:6)]
    fn <- unlist(strsplit(matrix_file, "/"))[length(unlist(strsplit(matrix_file, "/")))]
    runList <- list(runDef = runDef,
                    computeMatrix = computeMatrix,
                    computeMatrixRows = rows,
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



