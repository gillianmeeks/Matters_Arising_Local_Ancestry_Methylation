#' @export
flag.overlap <- function(probe_bed, SNP_bed) {
  # Recode chromosome information, if necessary
  if (length(grep("chr", probe_bed[,1], ignore.case=T)) > 0) {
    probe_bed[,1] <- gsub("chr", "", probe_bed[,1], ignore.case=T)
  }
  if (length(grep("chr", SNP_bed[,1], ignore.case=T)) > 0) {
    SNP_bed[,1] <- gsub("chr", "", SNP_bed[,1], ignore.case=T)
  }
  # Convert bed files to GRange objects
  subject <- GenomicRanges::GRanges(seqnames=probe_bed$chr, ranges=IRanges(start=probe_bed$start , end=probe_bed$end))
  #changed SNP range to only be one base as GenomicRanges uses closed intervals
  query <- GenomicRanges::GRanges(seqnames=SNP_bed[,1], ranges=IRanges(start=SNP_bed[,3], end=SNP_bed[,3]))
  # Find overlaps between SNPs and probes
  message("Calculating overlap between probe list and SNP list...")
  overlaps <- GenomicRanges::findOverlaps(query=query, subject=subject)
  probe_SNP_info <- as.data.frame(matrix(nrow=length(overlaps), ncol=12))
  colnames(probe_SNP_info) <- c("chr", "SNP_start", "SNP_end", "SNP_id", "SNP_ref", "SNP_alt", "CpG_start", "CpG_end", "CpG_id", "CpG_strand", "CpG_type", "CpG_pos")
  probe_SNP_info[,c("chr", "SNP_start", "SNP_end", "SNP_id", "SNP_ref", "SNP_alt")] <- SNP_bed[queryHits(overlaps),]
  probe_SNP_info[,c("CpG_start", "CpG_end", "CpG_id", "CpG_strand", "CpG_type", "CpG_pos")] <- probe_bed[subjectHits(overlaps), c("start", "end", "CpG_id", "strand", "type", "CpG_pos")]
  # Calculate SNP-CpG distance
  probe_SNP_info$SNP_CpG_distance <- probe_SNP_info$SNP_end - probe_SNP_info$CpG_pos
  # Annotate colour channel switching SNPs
  message("Annotating colour channel switching SNPs...")
  probe_SNP_info$col_chan_switching <- apply(X=probe_SNP_info, MARGIN=1, FUN=check.for.cc.switch)
  return(probe_SNP_info)
}