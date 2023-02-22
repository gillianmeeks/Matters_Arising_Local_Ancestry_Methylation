check.for.cc.switch <- function(DAT) {
  ref <- DAT[5]
  alt <- DAT[6]
  strand <- DAT[10]
  type <- DAT[11]
  distance <- as.numeric(DAT[13])
  if (type == 1) {
    if ((strand == '+' & distance == -1) | (strand == '-' & distance == 2)) {
      if ((ref == 'A' & alt == 'T') | (ref == 'T' & alt == 'A') | (ref == 'C' & alt == 'G') | (ref == 'G' & alt == 'C')) {
        return(F)
      } else {
        return(T)
      }
    } else {
      #added desginations for non-type 1 SBE SNPs
      return("not SBE")
    }
  } else {
    return("not SBE")
  }
}