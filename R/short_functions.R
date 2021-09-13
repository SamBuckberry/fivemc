

#' gr_to_loci
#'
#' @param gr GRanges object
#' @description Convert a GRanges object to a character vector of loci. These are useful as row.names for data.frame objects, pasting into a genome browser, among other things.
#' @return Character.
#' @export
#'
#' @examples gr <- GRanges(seqnames=Rle(c('chr1'), c(3)), IRanges(1:3, width=5))
#' gr_to_loci(gr)
#'
gr_to_loci <- function(gr){
    str_c(str_c(GenomicRanges::seqnames(gr), start(gr), sep = ":") %>%
              str_c(end(gr), sep = "-"))
}



