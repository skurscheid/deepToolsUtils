DownloadEnsemblGTF <- function(version = "current",
                                 url = "ftp.ensembl.org/pub/current_gtf",
                                 temp_dir = "/tmp",
                                 organism = "homo_sapiens"){
  match.arg(organism)
  source_file = paste(url, organism, "Homo_sapiens.GRCh38.84.gtf.gz", sep = "/")
  dest_file = paste(temp_dir, "Homo_sapiens.GRCh38.84.gtf.gz", sep = "/")
  if (! file.exists(dest_file)){
    tryCatch(curl::curl_download(url = source_file, destfile = dest_file, mode = "wb"))
    return(dest_file)
  } else {
    return(dest_file)
  }
}
