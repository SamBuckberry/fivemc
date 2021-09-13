

#' chDMR
#' @description Function for detecting CH (non-CG) differentially methylated regions (DMRs)
#' @param pData a `data.frame` object with the columns "sample", "group" and "bigwig_path". The "sample" column contains the unique identifiers for each sample. The "group" column contains the identifiers used to define groups for differential testing. The "bigwig_path" column contains the absolute file path for bigwig files created using the `bsseq_to_window_bigwig` function.
#' @param bsgenome
#' @param contigs
#' @param cores
#' @param normalize
#' @param normalize_method
#'
#' @return
#' @export
#'
#' @examples
chDMW <- function(pData, p=0.05,
                        bsgenome=BSgenome.Hsapiens.UCSC.hg19,
                        contigs=paste0("chr", c(1:22)), cores=2L,
                        normalize=TRUE, normalize_method="quantile"){

    # Check inputs
    assertthat::assert_that(is.data.frame(pData))
    assertthat::assert_that(all(c("sample", "group", "bigwig_path") %in%
                                    colnames(pData)))
    assertthat::assert_that(all(file.exists(pData$bigwig_path)))
    assertthat::assert_that(is.integer(cores))

    # Read the data from bigwig files
    read_bigwig <- function(bigwig_path, contig){

        # Get the info for the contigs that will be profiled
        gr_genome <- GenomicRanges::GRanges(seqinfo(BSgenome.Hsapiens.UCSC.hg19))
        gr_genome <- gr_genome[seqnames(gr_genome) == contig]

        bw_gr <- suppressWarnings(
            rtracklayer::import(con = bigwig_path,
                                format="bigWig",
                                selection = BigWigSelection(ranges = gr_genome),
                                as = c("GRanges"))
        )

        loci <- gr_to_loci(bw_gr)

        mat <- data.frame(loci, bw_gr$score)
        colnames(mat) <- c("loci", basename(bigwig_path))

        return(mat)
    }

    process_contig <- function(contig){

        message(paste0("Processing ", contig, ": ", as.character(Sys.time())))

        ch_dat <- parallel::mclapply(pData$bigwig_path,
                                     read_bigwig, contig=contig, mc.cores = cores)

        ch_dat <- base::Reduce(function(x, y) merge(x, y,
                                                    by = "loci", all = TRUE), ch_dat)
        row.names(ch_dat) <- ch_dat$loci
        ch_dat$loci <- NULL

        colnames(ch_dat) <- pData$sample

        return(ch_dat)
    }

    dat <- lapply(X = contigs, process_contig)
    dat <- do.call(rbind, dat)

    # Perform DMR testing

    groups <- factor(pData$group, levels=unique(pData$group))

    design <- model.matrix(~0+groups)
    colnames(design) <- levels(groups)

    # Eliminate rows with NA (indicative of no coverage)
    keep <- complete.cases(dat)
    dat <- dat[keep, ]

    # Eliminate windows with zero mCH
    keep2 <- rowSums(dat) > 0
    dat <- dat[keep2, ]

    # Normalise data
    if (normalize == TRUE) {
        message("Normalizing data using method: ", normalize_method, as.character(Sys.time()))
        dat <- limma::normalizeBetweenArrays(log2(dat+1),
                                      method = normalize_method)
    }

    message("Testing for differentially methylated windows: ", as.character(Sys.time()))
    fit <- limma::lmFit(dat, design)
    cont.matrix <- limma::makeContrasts(iPSCvsESC=iPSC-ESC, levels=design)
    fit2 <- limma::contrasts.fit(fit, cont.matrix)
    fit2 <- limma::eBayes(fit2)
    tt <- limma::topTable(fit2, adjust="BH", number = nrow(dat))

    dmw <- GRanges(rownames(tt))
    mcols(dmw) <- tt

    message(paste0("chDMW complete: ", as.character(Sys.time())))

    return(dmw)
}

#' chDMR
#' @description Find contiguous regions of significant changes with same direction
#' @param dmw A `GRanges` object output from the `chDMW` function.
#' @param p A `numeric` P-value threshold for selecting candidate differentially methylated windows. Default is 0.05.
#' @param min_size An `integer` of minimum CH-DMR size.
#' @param gap1_merge An `integer` representing the distance to merge deferentially methylated windows to create DMRs.
#' @param gap2_merge An `integer` representing the distance to merge candidate DMRs.
#' @return An object of class `GRanges`
#' @export
#'
#' @examples
chDMR <- function(dmw, p=0.05, min_size=60000,
                  gap1_merge=5000,
                  gap2_merge=100000){

    window_up <- dmw[dmw$P.Value < p & dmw$logFC > 0]
    window_up <- GenomicRanges::reduce(window_up, min.gapwidth=gap1_merge)
    window_up <- window_up[width(window_up) >= min_size]
    window_up_merged <- GenomicRanges::reduce(window_up, min.gapwidth=gap2_merge)

    window_down <- dmw[dmw$P.Value < p & dmw$logFC < 0]
    window_down <- GenomicRanges::reduce(window_down, min.gapwidth=gap1_merge)
    window_down <- window_down[width(window_down) >= min_size]
    window_down_merged <- GenomicRanges::reduce(window_down, min.gapwidth=gap2_merge)

    window_up_merged <- window_up_merged[!overlapsAny(window_up_merged, window_down)]
    window_down_merged <- window_down_merged[!overlapsAny(window_down_merged, window_up)]

    window_up_merged$direction <- "Up"
    window_down_merged$direction <- "Down"

    overlapsAny(window_up_merged, window_down_merged) %>% table()

    all_dmr <- c(window_up_merged, window_down_merged)

    return(all_dmr)
}




