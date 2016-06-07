# make_bed_from_ensembl_gtf.R

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
