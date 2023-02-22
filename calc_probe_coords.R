calc.probe.coords <- function(DAT) {
  strand <- DAT[4]
  type <- as.numeric(DAT[5])
  bp_pos <- as.numeric(DAT[3])
  if (strand == "+" & type == 1) {
    start <- bp_pos - 1
    end <- bp_pos + 49
  } else if (strand == "+" & type == 2) {
    start <- bp_pos
    end <- bp_pos + 50
  } else if (strand == "-" & type == 1) {
    start <- bp_pos - 48
    end <- bp_pos + 2
  } else if (strand == "-" & type == 2) {
    start <- bp_pos - 49
    end <- bp_pos + 1
  }
  #changed from open, close to close, close intervals
  return(c(start, end))
}