WriteGRangesToBED <- function(gr = NULL,
                                 out_file = NULL){
  tryCatch(if (is.null(gr)){
    print("GRanges object missing")
  } else if (is.null(out_file)) {
    print("Output file is missing")
  } else {
    df <- data.frame(seqnames = seqnames(gr),
                     starts = as(start(gr) - 1, "integer"),
                     ends = as(end(gr), "integer"),
                     names = names(gr),
                     scores = ".",
                     strand = strand(gr))
    write.table(df, file = out_file, quote = F, sep = "\t", row.names = F, col.names = F)
  })
}


