# Deterministic definition checks; no patient data and no clinical score inference.
a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
source(file.path(dirname(dirname(normalizePath(f))),"R/utils.R"))
root <- release_root()
m <- read_tsv(file.path(root,"resources/signatures/MES95.tsv"))
p <- read_tsv(file.path(root,"resources/signatures/Primary201_membership.tsv"))
stopifnot(nrow(m)==95L,!anyDuplicated(m$gene_symbol),nrow(p)==1729L,
          length(unique(p$model_instance_id))==201L)
cat("PASS: MES95 = 95 unique genes; Primary201 = 201 representations / 1729 membership rows.\n")

c<-read_tsv(file.path(root,"resources/signatures/CRC_BUD11.tsv"))
stopifnot(nrow(c)==11L,!anyDuplicated(c$canonical_gene))
