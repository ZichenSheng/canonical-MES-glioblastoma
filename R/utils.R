# Locate this repository without embedding machine paths.
release_root <- function() {
 a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=",a)])
 if (length(f)!=1L) stop("Run with Rscript")
 dirname(dirname(normalizePath(f)))
}
read_tsv <- function(path) read.delim(path, check.names=FALSE, stringsAsFactors=FALSE)
