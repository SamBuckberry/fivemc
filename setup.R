library(devtools)
library(usethis)
library(assertthat)

usethis::use_package("assertthat")
usethis::use_package("magrittr")
usethis::use_package("stringr")
usethis::use_package("BiocParallel")
usethis::use_package("BSgenome.Hsapiens.UCSC.hg19")
usethis::use_package("bsseq")
usethis::use_package("usethis")
usethis::use_package("rtracklayer")
usethis::use_package("GenomicRanges")
usethis::use_package("limma")

load_all()

document()
check()
install()

