#' @export
extract.probe.regions <- function(manifest_anno_object) {
  # Check input
  column_names <- colnames(manifest_anno_object)
  if (column_names[1] != "CpG_id") {
    stop("First column must be \"CpG_id\"")
  }
  if (column_names[2] != "chr") {
    stop("Second column must be \"chr\"")
  }
  if (column_names[3] != "CpG_pos") {
    stop("Third column must be \"CpG_pos\"")
  }
  if (column_names[4] != "strand") {
    stop("Fourth column must be \"strand\"")
  }
  if (column_names[5] != "type") {
    stop("Fifth column must be \"type\"")
  }

  # Recode chromosome information, if necessary
  if (length(grep("chr", manifest_anno_object$chr, ignore.case=T)) > 0) {
    manifest_anno_object$chr <- gsub("chr", "", manifest_anno_object$chr, ignore.case=T)
  }
  # Recode probe type information, if necessary
  if (length(grep("I", manifest_anno_object$type)) > 0) {
    manifest_anno_object$type[manifest_anno_object$type=="I"] <- 1
    manifest_anno_object$type[manifest_anno_object$type=="II"] <- 2
  }
  # Recode strand information, if necessary
  if (length(grep("F", manifest_anno_object$strand)) > 0) {
    manifest_anno_object$strand[manifest_anno_object$strand=="F"] <- '+'
  }
  if (length(grep("R", manifest_anno_object$strand)) > 0) {
    manifest_anno_object$strand[manifest_anno_object$strand=="R"] <- '-'
  }
  # Remove rows with no position information, if necessary
  if (length(which(is.na(manifest_anno_object$CpG_pos))) > 0) {
    manifest_anno_object <- manifest_anno_object[-which(is.na(manifest_anno_object$CpG_pos)),]
  }

  # Add probe coordinates
  manifest_anno_object[,c("start", "end")] <- t(apply(X=manifest_anno_object, MARGIN=1, FUN=calc.probe.coords))
  manifest_anno_object <- manifest_anno_object[,c("chr","start","end","CpG_id","strand","type","CpG_pos")]
  return(manifest_anno_object)
}
