CheckAvailableEnsemblReleases <- function(url = "ftp://ftp.ensembl.org/pub/"){
  releases <- getURL(url = url, dirlistonly = TRUE)
  releases <- strsplit(releases, "\r*\n")[[1]]
  releases <- releases[grep("release", releases)]
  return(releases)
}

CheckAvailableEnsemblOrganisms <- function(url = "ftp://ftp.ensembl.org/pub/", release = NULL){
  if (is.null(release)) stop("Ensembl release number missing. Please check with CheckAvailableEnsemblReleases() and specify in function call.")
  tryCatch(match.arg(release, choices = CheckAvailableEnsemblReleases()))
  organisms <- getURL(url = paste(url, release, "/fasta/", sep = ""), dirlistonly = TRUE)
  organisms <- strsplit(organisms, "\r*\n")[[1]]
  return(organisms)
}

DownloadEnsemblGTF <- function(version = "current",
                            url = "ftp.ensembl.org/pub/current_gtf",
                            temp_dir = "/tmp",
                            organism = "homo_sapiens",
                            data_type = "genes"){
  match.arg(organism, choices = CheckAvailableEnsemblOrganisms())
  match.arg(data_type, choices = c())

  source_file = paste(url, organism, "Homo_sapiens.GRCh38.84.gtf.gz", sep = "/")
  dest_file = paste(temp_dir, "Homo_sapiens.GRCh38.84.gtf.gz", sep = "/")

  if (! file.exists(dest_file)){
    tryCatch(curl::curl_download(url = source_file, destfile = dest_file, mode = "wb"))
    return(dest_file)
  } else {
    return(dest_file)
  }
}

MakeBEDFromEnsemblGTF <- function(gtf_file = NULL,
                                  output_dir = "/tmp",
                                  features = c("genes", "transcripts", "exons", "cds", "promoters")){
  match.arg(features, several.ok = FALSE)

  if (tryCatch(file.exists(gtf_file))){
    tryCatch(txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file))
    if (features == "genes"){
      gr <- genes(txdb)
    } else if (features == "transcripts") {
      gr <- transcripts(txdb)
    } else if (features == "exons") {
      gr <- exons(txdb)
    } else if (features == "cds") {
      gr <- cds(txdb)
    } else if (features == "promoters") {
      gr <- promoters(txdb)
    }
    gr <- sort(gr)

  } else {
    print("Please call download_ensembl_gtf() first, or provide full path to your GTF file!")
  }
}

WriteGRangesToBED <- function(gr = NULL,
                              out_file = NULL){
  if (is.null(gr)) {stop("GRanges object missing")}
  if (is.null(out_file)) {stop("Output file is missing")}

  df <- data.frame(seqnames = seqnames(gr),
                   starts = as(start(gr) - 1, "integer"),
                   ends = as(end(gr), "integer"),
                   names = names(gr),
                   scores = ".",
                   strand = strand(gr))
  write.table(df, file = out_file, quote = F, sep = "\t", row.names = F, col.names = F)
}


