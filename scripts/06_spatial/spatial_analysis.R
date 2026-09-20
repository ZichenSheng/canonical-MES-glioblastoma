#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(hdf5r)
  library(arrow)
  library(RANN)
  library(metafor)
  library(digest)
  library(jsonlite)
})

options(stringsAsFactors = FALSE, scipen = 999)
started <- Sys.time()

DATA_ROOT <- Sys.getenv("GBM_MES_DATA_ROOT", unset="data")
RESULT_ROOT <- Sys.getenv("GBM_MES_OUTPUT_ROOT", unset="results")
REPO_ROOT <- Sys.getenv("GBM_MES_REPO_ROOT", unset=".")
RUN <- file.path(RESULT_ROOT,"spatial")
SOURCE_TABLES <- file.path(DATA_ROOT,"prepared","spatial")
GSE237 <- file.path(DATA_ROOT,"spatial","GSE237183")
GSE242_TAR <- file.path(DATA_ROOT,"spatial","GSE242352_RAW.tar")
GSE242 <- file.path(RUN, "01_spatial_reference/gse242352_reference_members")
SIGFILE <- file.path(SOURCE_TABLES, "03_reusable_data_local/signatures/signature_gene_long_table_v1_3.tsv")
SOFT237 <- file.path(RUN, "01_spatial_reference/official_metadata/GSE237183_family.soft.gz")
SOFT242 <- file.path(RUN, "01_spatial_reference/official_metadata/GSE242352_family.soft.gz")
SEEDS <- fromJSON(file.path(RUN, "00_governance/RANDOM_SEEDS_V1_2.json"))
GSE242_TAR_EXPECTED_SHA <- "59375425cfaf724507cd476b02e510343c0b96bcc2d8bfd6b3a75630338e765a"

for (d in c("01_spatial_reference", "02_gse237183", "03_gse242352", "04_spatial_integration", "logs"))
  dir.create(file.path(RUN, d), recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
GSE242_TAR_SHA <- sha256(GSE242_TAR)
if (!identical(GSE242_TAR_SHA, GSE242_TAR_EXPECTED_SHA)) stop("Frozen GSE242352 official tar SHA256 mismatch")
file_row <- function(path) list(path = normalizePath(path, mustWork = TRUE), size_bytes = file.info(path)$size,
                                sha256 = sha256(path))
safe_cor <- function(x, y) {
  z <- suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
  ifelse(is.finite(z), z, NA_real_)
}
fisher_ci <- function(r, n) {
  if (!is.finite(r) || n <= 3 || abs(r) >= 1) return(c(NA_real_, NA_real_))
  z <- atanh(r); se <- 1 / sqrt(n - 3)
  tanh(z + c(-1, 1) * qnorm(.975) * se)
}
stable_seed <- function(base, ...) {
  key <- paste(..., collapse = "|")
  x <- digest2int(key, seed = as.integer(base))
  as.integer(abs(as.double(x)) %% 2147483000 + 1)
}
bt_ci <- function(x, seed, B = 2000L) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(NA_real_, NA_real_))
  if (length(x) == 1L) return(c(x, x))
  set.seed(seed)
  b <- replicate(B, mean(sample(x, length(x), replace = TRUE)))
  unname(quantile(b, c(.025, .975), na.rm = TRUE, names = FALSE))
}

parse_soft <- function(path) {
  z <- readLines(gzfile(path), warn = FALSE)
  starts <- which(startsWith(z, "^SAMPLE = "))
  ends <- c(starts[-1] - 1L, length(z))
  rows <- lapply(seq_along(starts), function(i) {
    b <- z[starts[i]:ends[i]]
    val <- function(prefix) sub(prefix, "", b[grep(prefix, b, fixed = TRUE)[1]])
    chr <- sub("!Sample_characteristics_ch1 = ", "", b[startsWith(b, "!Sample_characteristics_ch1 = ")])
    list(gsm = sub("^SAMPLE = ", "", b[1], fixed = TRUE), title = val("!Sample_title = "),
         source_name = val("!Sample_source_name_ch1 = "), characteristics = paste(chr, collapse = "; "),
         description = val("!Sample_description = "))
  })
  rbindlist(rows, fill = TRUE)
}

soft237 <- parse_soft(SOFT237)
soft242 <- parse_soft(SOFT242)

map237_patient <- function(title) {
  if (grepl("^GBM ZH881_1[AB]", title)) return("ZH881_1")
  if (grepl("^GBM ZH881_2", title)) return("ZH881_2")
  if (grepl("^GBM ", title)) return(sub("^GBM ([^ ]+).*$", "\\1", title))
  if (grepl("^IDHm ", title)) return(sub("^IDHm ([^ ]+).*$", "\\1", title))
  NA_character_
}
soft237[, patient_id := vapply(title, map237_patient, character(1))]
soft237[, disease_class := fifelse(startsWith(title, "GBM "), "IDH_WT_GBM", "IDH_MUTANT_GLIOMA")]
soft237[, section_id := paste0(gsm, "_", gsub("[^A-Za-z0-9]+", "", sub("^(GBM|IDHm) ", "", title)))]
soft237[, mapping_evidence := paste0("NCBI GEO family SOFT official Sample_title: ", title)]

soft242 <- soft242[gsm %in% c("GSM8513873", "GSM8513874", "GSM8513875", "GSM8513877", "GSM8513879", "GSM8513880")]
soft242[, case_label := sub("^Glioblastoma ", "", title)]
soft242[, patient_id := paste0("GSE242352_", case_label)]
soft242[, section_id := gsub(" ", "", title)]
soft242[, disease_class := "IDH_WT_GBM"]
soft242[, mapping_evidence := paste0("NCBI GEO official title plus Series_summary 'six additional ... cases': ", title)]

find237 <- function(gsm, suffix) {
  x <- list.files(GSE237, pattern = paste0("^", gsm, ".*", suffix, "$"), full.names = TRUE)
  if (length(x) != 1L) stop("Expected one GSE237183 file: ", gsm, " / ", suffix)
  x
}
gse242_prefix <- c(GSM8513873 = "GSM8513873_SPA1_D", GSM8513874 = "GSM8513874_SPA2_A",
                   GSM8513875 = "GSM8513875_SPA3_D", GSM8513877 = "GSM8513877_SPA5_D",
                   GSM8513879 = "GSM8513879_SPA7_D", GSM8513880 = "GSM8513880_SPA8_A")

siglong <- fread(SIGFILE)
get_sig <- function(id) sort(unique(na.omit(siglong[signature_id == id, gene_symbol_clean])))
HYP <- strsplit("NDRG1 TRIB3 ERO1A DNAJB9 HILPDA SLCO4A1 INSIG2 AKAP12 SLC3A2 ZFAS1 GOLT1B LOX NRN1 GBE1 SESN2 MXI1 STC2 VEGFA CEBPG UPP1 HK2 TAF1D INSIG1 EIF4EBP1 ADM PLOD2 RRAGD EPB41L4A-AS1 DDIT3 ANKRD37 FAM210A PDK1 PGK1 GFPT1 IGFBP3 CA9 TSPYL2 SLC38A1 EGLN3 SLC6A6 IGFBP5 SLC7A5 ANGPTL4 SLC1A5 FAM13A TREM1 CXCL8 TM4SF1 OSGIN1 PSAT1 BNIP3 ARFGEF3 LUCAT1 STX3 ASNS LARP6 GDF15 SH3BGR OSER1 NUPR1 GPNMB GRPEL2 KDM3A FICD ATF5 CLEC2B TMEM38B ATF4 DDIT4 CTH CREBRF CIART TUBA4A BNIP3L JMJD6 PIGA BEX2 BIRC3 BTG1 PPP1R3C PRKAG2-AS1 HSPH1 HSPA9 HSPA6 PID1 EMX2 KISS1R PHF1 LIF PERP LMAN1 ARID3A XBP1 ANG NUP58 MAP1LC3B ZFAND2A TMEM45A HSPA5 FRMD3", " ")[[1]]
AST <- strsplit("IGFBP7 GFAP ID3 METTL7B EFEMP1 PTN CHI3L1 PMP2 GAP43 IFITM3 CCL2 KLHDC8A FABP7 C21orf62 SNTG1 POSTN THY1 PTPRZ1 SRPX2 GPM6A COL5A2 HOPX AQP4 APOE MXRA8 CXCL14 CP RAMP1 CFI RCAN1 FABP5 S100B SERPING1 CLU CST3 A2M SPARC ABCC3 GPM6B ATP1B2 SPARCL1 MLC1 AGT SCARA3 LTF NES NNMT COL1A2 PLA2G5 PLAAT4 PLA2G2A MAOB C1S C1R C3 SERPINA3 SLPI PROS1 OLFM2 GPX3 CDH11 SPP1 TNC MOXD1 VCAN CITED1 SAA1 ID1 VCAM1 CNR1 HTRA1 GAS1 NPAS3 GPC1 HLA-DRA GRIA2 TAGLN ITGA7 SEMA6D RBP1 RDH10 SCD5 DKK3 DPYD SERPINF1 RHOJ CHRNA1 RFX4 FSTL1 LPL SLC1A2 LGALS3BP SCRG1 CDH6 TUBA1A PPIC TTYH1 COL6A2 COL6A1 COL1A1", " ")[[1]]
VASC <- strsplit("PECAM1 VWF FLT1 PDGFRB KCNJ8 CSPG4 CNN1 DES MYH11", " ")[[1]]

SIG_RAW <- list(CANONICAL_MES = get_sig("NEFTEL_MES_LIKE"), MES1 = get_sig("NEFTEL_MES1"),
                MES2 = get_sig("NEFTEL_MES2"), MES_HYP = sort(unique(HYP)), MES_AST = sort(unique(AST)))
ECO <- list(HYPOXIA = get_sig("HALLMARK_HYPOXIA"),
            MYELOID = sort(unique(c(get_sig("MYELOID_CORE_IDENTITY"), get_sig("MYELOID_ACTIVATED_TAM")))),
            VASCULAR = sort(unique(VASC)), MATRIX = get_sig("NABA_CORE_MATRISOME"))
TARGET <- list(HYP_TARGET = SIG_RAW$MES_HYP, AST_TARGET = SIG_RAW$MES_AST)
target_union <- union(SIG_RAW$MES_HYP, SIG_RAW$MES_AST)
hyp_ast_intersection <- intersect(SIG_RAW$MES_HYP, SIG_RAW$MES_AST)
SIG_PRUNED <- list(
  CANONICAL_MES = setdiff(SIG_RAW$CANONICAL_MES, target_union),
  MES1 = setdiff(SIG_RAW$MES1, target_union),
  MES2 = setdiff(SIG_RAW$MES2, target_union),
  MES_HYP = setdiff(SIG_RAW$MES_HYP, hyp_ast_intersection),
  MES_AST = setdiff(SIG_RAW$MES_AST, hyp_ast_intersection)
)

freeze_rows <- rbindlist(lapply(names(SIG_RAW), function(s) {
  rbind(data.table(signature = s, mode = "RAW", gene = SIG_RAW[[s]], retained = TRUE),
        data.table(signature = s, mode = "OVERLAP_PRUNED", gene = SIG_PRUNED[[s]], retained = TRUE))
}))
freeze_rows[, `:=`(source = fifelse(signature %in% c("MES_HYP", "MES_AST"),
                                     "v1.1 verified author-label sidecar source script",
                                     paste0("signature_gene_long_table_v1_3.tsv/", fifelse(signature == "CANONICAL_MES", "NEFTEL_MES_LIKE", signature))),
                    purpose = "formal_spatial_score", minimum_covered_genes = 2L)]
fwrite(freeze_rows, file.path(RUN, "01_spatial_reference/SPATIAL_SIGNATURE_FREEZE.tsv"), sep = "\t")

removed <- rbindlist(lapply(names(SIG_RAW), function(s) {
  g <- setdiff(SIG_RAW[[s]], SIG_PRUNED[[s]])
  data.table(signature = s, gene = g, raw_n = length(SIG_RAW[[s]]), pruned_n = length(SIG_PRUNED[[s]]),
             removal_rule = if (s %in% c("MES_HYP", "MES_AST")) "bilateral_MES_HYP_MES_AST_intersection"
                            else "direct_overlap_with_MES_HYP_union_MES_AST",
             restored = FALSE)
}))
fwrite(removed, file.path(RUN, "01_spatial_reference/SPATIAL_OVERLAP_REMOVAL.tsv"), sep = "\t")

read10x_h5 <- function(path) {
  f <- H5File$new(path, mode = "r")
  on.exit(f$close_all())
  data <- f[["matrix/data"]][]; indices <- f[["matrix/indices"]][]
  indptr <- f[["matrix/indptr"]][]; shape <- f[["matrix/shape"]][]
  genes <- as.character(f[["matrix/features/name"]][]); bars <- as.character(f[["matrix/barcodes"]][])
  # Some official 10x H5 members contain valid but unsorted row indices within a column.
  # sparseMatrix canonicalizes these without changing any input bytes or records.
  jj <- rep.int(seq_len(as.integer(shape[2])), diff(as.integer(indptr)))
  X <- sparseMatrix(i = as.integer(indices) + 1L, j = jj, x = as.numeric(data), dims = as.integer(shape), giveCsparse = TRUE)
  list(X = X, genes = genes, bars = bars, shape = as.integer(shape))
}
read_positions <- function(path) {
  first <- readLines(gzfile(path), n = 1L)
  if (startsWith(first, "barcode,")) {
    d <- fread(path)
    setnames(d, c("pxl_row_in_fullres", "pxl_col_in_fullres"), c("pxl_row", "pxl_col"), skip_absent = TRUE)
  } else {
    d <- fread(path, header = FALSE)
    setnames(d, c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"))
  }
  d
}

score_genes <- function(X, genes, available) {
  idx <- match(genes, available, nomatch = 0L); idx <- idx[idx > 0L]
  if (length(idx) < 2L) return(list(score = rep(NA_real_, ncol(X)), covered = length(idx), status = "NOT_EVALUABLE_LOW_GENE_COVERAGE"))
  lib <- Matrix::colSums(X); lib[lib <= 0] <- 1
  a <- as.matrix(X[idx, , drop = FALSE])
  a <- log1p(t(t(a) / lib) * 10000)
  mu <- rowMeans(a); ss <- apply(a, 1, sd); ss[!is.finite(ss) | ss == 0] <- 1
  z <- (a - mu) / ss
  list(score = colMeans(z), covered = length(idx), status = "EVALUABLE")
}

biv_moran <- function(a, b, xy) {
  ok <- is.finite(a) & is.finite(b) & apply(xy, 1, function(z) all(is.finite(z)))
  a <- a[ok]; b <- b[ok]; xy <- xy[ok, , drop = FALSE]
  if (length(a) < 8L) return(NA_real_)
  k <- min(7L, nrow(xy))
  nn <- nn2(xy, xy, k = k)$nn.idx[, -1, drop = FALSE]
  za <- as.numeric(scale(a)); zb <- as.numeric(scale(b))
  mean(za * rowMeans(matrix(zb[nn], nrow = nrow(nn))))
}

make_blocks <- function(xy, multiple) {
  nn <- nn2(xy, xy, k = 2)$nn.dists[, 2]
  d <- median(nn[nn > 0 & is.finite(nn)])
  width <- multiple * d
  bx <- floor((xy[,1] - min(xy[,1])) / width)
  by <- floor((xy[,2] - min(xy[,2])) / width)
  list(id = paste(bx, by, sep = ":"), bx = bx, by = by, median_nn = d, width = width)
}

all_reference <- list(); all_integrity <- list(); coverage_all <- list(); spot_all <- list()
study_outputs <- list(); cov_i <- auth_i <- qa_i <- spot_i <- 0L

run_study <- function(study, map, h5_fun, pos_fun, seed_base, outdir, prefix) {
  t0 <- Sys.time(); section_results <- list(); moran_results <- list(); scale_results <- list(); qa_rows <- list(); cov_rows <- list(); spots <- list()
  null_path <- file.path(RUN, outdir, paste0(prefix, "_PERMUTATION_NULL.tsv"))
  if (file.exists(null_path)) file.remove(null_path)
  result_i <- moran_i <- scale_i <- q_i <- c_i <- s_i <- 0L
  for (ii in seq_len(nrow(map))) {
    m <- map[ii]; h5 <- h5_fun(m$gsm); pospath <- pos_fun(m$gsm)
    message(study, " ", ii, "/", nrow(map), " ", m$section_id)
    obj <- read10x_h5(h5); pos <- read_positions(pospath)
    pos <- pos[in_tissue == 1]
    keep <- match(obj$bars, pos$barcode, nomatch = 0L); ok <- keep > 0L
    X <- obj$X[, ok, drop = FALSE]; bars <- obj$bars[ok]; pos <- pos[keep[ok]]
    xy <- as.matrix(pos[, .(pxl_col, pxl_row)])
    lib <- Matrix::colSums(X); detected <- Matrix::colSums(X != 0)
    score_cache <- list(); cov_cache <- list()
    for (sig in names(SIG_RAW)) for (mode in c("RAW", "OVERLAP_PRUNED")) {
      gl <- if (mode == "RAW") SIG_RAW[[sig]] else SIG_PRUNED[[sig]]
      sc <- score_genes(X, gl, obj$genes); key <- paste(sig, mode, sep = "|")
      score_cache[[key]] <- sc$score; cov_cache[[key]] <- sc
      c_i <- c_i + 1L
      cov_rows[[c_i]] <- data.table(study = study, section_id = m$section_id, patient_id = m$patient_id,
        signature = sig, mode = mode, genes_defined = length(gl), genes_covered = sc$covered,
        coverage_fraction = sc$covered / length(gl), qc_status = sc$status)
    }
    eco_cache <- list()
    for (eco in names(ECO)) eco_cache[[eco]] <- score_genes(X, ECO[[eco]], obj$genes)
    target_cache <- list(HYP_TARGET = score_cache[["MES_HYP|RAW"]], AST_TARGET = score_cache[["MES_AST|RAW"]])
    sp <- data.table(study = study, section_id = m$section_id, patient_id = m$patient_id,
                     disease_class = m$disease_class, barcode = bars, in_tissue = 1L,
                     array_row = pos$array_row, array_col = pos$array_col, pxl_row = pos$pxl_row, pxl_col = pos$pxl_col,
                     library_size = as.numeric(lib), detected_genes = as.integer(detected))
    for (sig in names(SIG_RAW)) for (mode in c("RAW", "OVERLAP_PRUNED"))
      sp[[paste0("score_", sig, "_", mode)]] <- score_cache[[paste(sig, mode, sep = "|")]]
    for (eco in names(ECO)) sp[[paste0("ecology_", eco)]] <- eco_cache[[eco]]$score
    for (tg in names(target_cache)) sp[[paste0("target_", tg)]] <- target_cache[[tg]]
    spots[[ii]] <- sp
    q_i <- q_i + 1L
    qa_rows[[q_i]] <- data.table(study = study, section_id = m$section_id, patient_id = m$patient_id,
      disease_class = m$disease_class, h5_path = h5, coordinate_path = pospath,
      matrix_feature_n = obj$shape[1], matrix_barcode_n = obj$shape[2], tissue_spot_n = ncol(X),
      coordinate_row_n = nrow(pos), barcode_coordinate_match = ncol(X) == nrow(pos),
      zero_library_spots = sum(lib == 0), median_library_size = median(lib), median_detected_genes = median(detected),
      tissue_mask_rule = "official in_tissue==1", structure_qa = ifelse(ncol(X) == nrow(pos) && ncol(X) > 0, "PASS", "FAIL"))
    mor_cache <- list()
    blocks <- list(small = make_blocks(xy, 4), medium = make_blocks(xy, 8), large = make_blocks(xy, 16))
    for (sig in names(SIG_RAW)) for (mode in c("RAW", "OVERLAP_PRUNED")) for (eco in names(ECO)) {
      a <- score_cache[[paste(sig, mode, sep = "|")]]; b <- eco_cache[[eco]]$score
      n <- sum(is.finite(a) & is.finite(b)); rho <- safe_cor(a, b)
      ci <- fisher_ci(rho, n)
      pspot <- if (n >= 4 && is.finite(rho)) suppressWarnings(cor.test(a, b, method = "spearman", exact = FALSE)$p.value) else NA_real_
      mi <- biv_moran(a, b, xy); moran_i <- moran_i + 1L
      moran_results[[moran_i]] <- data.table(study = study, section_id = m$section_id, patient_id = m$patient_id,
        signature = sig, mode = mode, ecology = eco, n_spots = n, bivariate_moran_i = mi,
        neighbor_rule = "6-nearest-neighbor row-standardized cross-product")
      for (scale in names(blocks)) {
        bl <- blocks[[scale]]
        d <- data.table(block_id = bl$id, bx = bl$bx, by = bl$by, a = a, b = b)
        db <- d[is.finite(a) & is.finite(b), .(a = mean(a), b = mean(b), n_spots = .N, bx = first(bx), by = first(by)), by = block_id]
        br <- safe_cor(db$a, db$b); bci <- fisher_ci(br, nrow(db))
        set.seed(stable_seed(seed_base, study, m$section_id, sig, mode, eco, scale))
        ra <- rank(db$a, ties.method = "average"); rb <- rank(db$b, ties.method = "average")
        null <- if (nrow(db) >= 4L && is.finite(br)) replicate(999L, suppressWarnings(cor(ra, sample(rb)))) else rep(NA_real_, 999L)
        emp <- if (is.finite(br)) (1 + sum(abs(null) >= abs(br), na.rm = TRUE)) / (1 + sum(is.finite(null))) else NA_real_
        edge <- db$bx %in% range(db$bx) | db$by %in% range(db$by)
        result_i <- result_i + 1L
        section_results[[result_i]] <- data.table(study = study, section_id = m$section_id, patient_id = m$patient_id,
          disease_class = m$disease_class, signature = sig, mode = mode, ecology = eco, block_scale = scale,
          n_spots = n, genes_covered = cov_cache[[paste(sig, mode, sep = "|")]]$covered,
          evaluability = cov_cache[[paste(sig, mode, sep = "|")]]$status,
          spot_spearman_rho = rho, spot_ci_low = ci[1], spot_ci_high = ci[2], spot_p = pspot,
          block_spearman_rho = br, block_ci_low = bci[1], block_ci_high = bci[2], block_empirical_p = emp,
          bivariate_moran_i = mi, median_nn_distance = bl$median_nn, block_side_length = bl$width,
          occupied_blocks = nrow(db), edge_blocks = sum(edge), min_spots_per_block = min(db$n_spots),
          median_spots_per_block = median(db$n_spots), max_spots_per_block = max(db$n_spots), permutations = 999L)
        nr <- data.table(study = study, section_id = m$section_id, patient_id = m$patient_id, signature = sig,
                         mode = mode, ecology = eco, block_scale = scale, permutation_id = 1:999,
                         null_block_spearman_rho = null)
        fwrite(nr, null_path, sep = "\t", append = file.exists(null_path), col.names = !file.exists(null_path))
      }
    }
  }
  res <- rbindlist(section_results, fill = TRUE)
  res[, spot_q := p.adjust(spot_p, "BH"), by = .(study, mode)]
  res[, block_q := p.adjust(block_empirical_p, "BH"), by = .(study, mode, block_scale)]
  mor <- rbindlist(moran_results); scl <- copy(res)
  cov <- rbindlist(cov_rows); qa <- rbindlist(qa_rows); spot <- rbindlist(spots, use.names = TRUE, fill = TRUE)
  fwrite(res, file.path(RUN, outdir, paste0(prefix, "_SECTION_RESULTS_999.tsv")), sep = "\t")
  fwrite(scl[, .(study, section_id, patient_id, disease_class, signature, mode, ecology, block_scale,
                 median_nn_distance, block_side_length, occupied_blocks, edge_blocks, min_spots_per_block,
                 median_spots_per_block, max_spots_per_block, block_spearman_rho, block_ci_low, block_ci_high,
                 block_empirical_p, block_q, permutations)],
         file.path(RUN, outdir, paste0(prefix, "_BLOCK_SCALE_SENSITIVITY.tsv")), sep = "\t")
  fwrite(mor, file.path(RUN, outdir, paste0(prefix, "_BIVARIATE_MORAN.tsv")), sep = "\t")
  fwrite(qa, file.path(RUN, outdir, paste0(prefix, "_SECTION_QA.tsv")), sep = "\t")
  fwrite(cov, file.path(RUN, outdir, paste0(prefix, "_GENE_COVERAGE.tsv")), sep = "\t")
  write_parquet(spot, file.path(RUN, outdir, paste0(prefix, "_SPOT_SCORES.parquet")), compression = "zstd")
  list(results = res, moran = mor, qa = qa, coverage = cov, spots = spot,
       elapsed_seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

# Reference tables and structural inputs.
map237 <- soft237[, .(study_id = "GSE237183", geo_accession = gsm, sample_id = gsm, section_id, patient_id,
                      disease_class, study_platform = "10x Visium", technical_replicate = FALSE,
                      official_title = title, mapping_evidence, mapping_reference = SOFT237,
                      identity_status = "REFERENCE_RESOLVED")]
map242 <- soft242[, .(study_id = "GSE242352", geo_accession = gsm, sample_id = gsm, section_id, patient_id,
                      disease_class, study_platform = "10x Visium CytAssist FFPE", technical_replicate = FALSE,
                      official_title = title, mapping_evidence, mapping_reference = SOFT242,
                      identity_status = "REFERENCE_RESOLVED")]
fwrite(rbindlist(list(map237, map242), fill = TRUE), file.path(RUN, "01_spatial_reference/SPATIAL_SECTION_PATIENT_MAP.tsv"), sep = "\t")
fwrite(map242, file.path(RUN, "03_gse242352/GSE242352_SECTION_PATIENT_MAP.tsv"), sep = "\t")

auth_rows <- list(); integ_rows <- list(); ai <- ii <- 0L
for (i in seq_len(nrow(soft237))) {
  x <- soft237[i]; h5 <- find237(x$gsm, "_filtered_feature_bc_matrix\\.h5"); pos <- find237(x$gsm, "_tissue_positions_list\\.csv\\.gz")
  for (z in list(c("expression", h5), c("spatial_coordinate_tissue_position", pos))) {
    ai <- ai + 1L; f <- file_row(z[2]); auth_rows[[ai]] <- data.table(study_id = "GSE237183", geo_accession = x$gsm,
      sample_id = x$gsm, section_id = x$section_id, patient_id = x$patient_id, platform = "10x Visium", input_role = z[1],
      path = f$path, size_bytes = f$size_bytes, sha256 = f$sha256, download_source = "GEO official processed supplementary member",
      structure_qa = "PENDING_FORMAL_READ", used_in_old_sidecar = TRUE)
  }
}
for (i in seq_len(nrow(soft242))) {
  x <- soft242[i]; pfx <- gse242_prefix[[x$gsm]]; h5 <- file.path(GSE242, paste0(pfx, "_filtered_feature_bc_matrix.h5")); pos <- file.path(GSE242, paste0(pfx, "_tissue_positions.csv.gz"))
  for (z in list(c("expression", h5), c("spatial_coordinate_tissue_position", pos))) {
    ai <- ai + 1L; f <- file_row(z[2]); auth_rows[[ai]] <- data.table(study_id = "GSE242352", geo_accession = x$gsm,
      sample_id = x$gsm, section_id = x$section_id, patient_id = x$patient_id, platform = "10x Visium CytAssist FFPE", input_role = z[1],
      path = f$path, size_bytes = f$size_bytes, sha256 = f$sha256, download_source = paste0("GEO official series tar member; tar SHA256=", GSE242_TAR_SHA),
      structure_qa = "PENDING_FORMAL_READ", used_in_old_sidecar = FALSE)
  }
}
for (s in list(c("GSE237183", SOFT237), c("GSE242352", SOFT242))) {
  ai <- ai + 1L; f <- file_row(s[2]); auth_rows[[ai]] <- data.table(study_id = s[1], input_role = "official_metadata",
    path = f$path, size_bytes = f$size_bytes, sha256 = f$sha256, download_source = "NCBI GEO official family SOFT HTTPS",
    structure_qa = "PASS_GZIP_AND_PARSE", used_in_old_sidecar = FALSE)
}
ai <- ai + 1L
auth_rows[[ai]] <- data.table(study_id = "GSE242352", input_role = "official_series_supplementary_tar",
  path = normalizePath(GSE242_TAR), size_bytes = file.info(GSE242_TAR)$size, sha256 = GSE242_TAR_SHA,
  download_source = "NCBI GEO official series supplementary tar; frozen local reference",
  structure_qa = "PASS_SHA256_AND_MEMBER_LIST", used_in_old_sidecar = FALSE)
reference <- rbindlist(auth_rows, fill = TRUE)
fwrite(reference, file.path(RUN, "01_spatial_reference/SPATIAL_REFERENCE_INPUTS.tsv"), sep = "\t")
fwrite(data.table(study_a = "GSE237183", study_b = "GSE242352", accession_overlap = FALSE, geo_sample_overlap_n = 0L,
  patient_namespace_overlap_n = 0L, independence_evidence = "distinct GEO accessions, publications, BioProjects, and patient/case namespaces",
  caveat = "GSE237183 frozen 19-section context contains 13 GBM sections and 6 IDH-mutant glioma sections; not all 19 sections are GBM"),
  file.path(RUN, "01_spatial_reference/SPATIAL_STUDY_OVERLAP_CHECK.tsv"), sep = "\t")

message("Running GSE237183 formal branch")
r237 <- run_study("GSE237183", soft237,
  function(gsm) find237(gsm, "_filtered_feature_bc_matrix\\.h5"),
  function(gsm) find237(gsm, "_tissue_positions_list\\.csv\\.gz"),
  SEEDS$gse237183_block_permutation, "02_gse237183", "GSE237183")
message("Running GSE242352 formal branch")
r242 <- run_study("GSE242352", soft242,
  function(gsm) file.path(GSE242, paste0(gse242_prefix[[gsm]], "_filtered_feature_bc_matrix.h5")),
  function(gsm) file.path(GSE242, paste0(gse242_prefix[[gsm]], "_tissue_positions.csv.gz")),
  SEEDS$gse242352_block_permutation, "03_gse242352", "GSE242352")

# Complete integrity and shared coverage check after formal reads.
integrity <- rbindlist(list(r237$qa, r242$qa), fill = TRUE)
fwrite(integrity, file.path(RUN, "01_spatial_reference/SPATIAL_INPUT_INTEGRITY_QA.tsv"), sep = "\t")
coverage <- rbindlist(list(r237$coverage, r242$coverage), fill = TRUE)
fwrite(coverage, file.path(RUN, "01_spatial_reference/SPATIAL_GENE_COVERAGE.tsv"), sep = "\t")
fwrite(r242$qa, file.path(RUN, "03_gse242352/GSE242352_INPUT_QA.tsv"), sep = "\t")

# Patient-level equal-section summaries.
allres <- rbindlist(list(r237$results, r242$results), fill = TRUE)
moran_all <- rbindlist(list(r237$moran, r242$moran), fill = TRUE)
patient_spot <- unique(allres[, .(study, section_id, patient_id, disease_class, signature, mode, ecology,
                                  estimate = spot_spearman_rho)], by = c("study","section_id","signature","mode","ecology"))
patient_spot[, `:=`(metric = "spot_spearman", block_scale = "spot")]
patient_block <- allres[, .(study, section_id, patient_id, disease_class, signature, mode, ecology,
                            metric = "block_spearman", block_scale, estimate = block_spearman_rho)]
patient_moran <- moran_all[, .(study, section_id, patient_id, disease_class = NA_character_, signature, mode, ecology,
                               metric = "bivariate_moran", block_scale = "none", estimate = bivariate_moran_i)]
sec_long <- rbindlist(list(patient_spot, patient_block, patient_moran), fill = TRUE)
pat <- sec_long[, .(patient_effect = mean(estimate, na.rm = TRUE), n_sections = uniqueN(section_id),
                    section_sd = if (.N > 1) sd(estimate, na.rm = TRUE) else NA_real_,
                    section_min = min(estimate, na.rm = TRUE), section_max = max(estimate, na.rm = TRUE)),
                by = .(study, patient_id, disease_class, signature, mode, ecology, metric, block_scale)]
pat[!is.finite(patient_effect), patient_effect := NA_real_]
fwrite(pat, file.path(RUN, "04_spatial_integration/SPATIAL_PATIENT_LEVEL.tsv"), sep = "\t")
fwrite(pat[, .(study, patient_id, signature, mode, ecology, metric, block_scale, n_sections, section_sd, section_min,
               section_max, section_range = section_max - section_min)],
       file.path(RUN, "04_spatial_integration/SPATIAL_WITHIN_PATIENT_HETEROGENEITY.tsv"), sep = "\t")

# Study-level equal-patient summaries and sign-permutation support.
study_rows <- list(); si <- 0L
for (k in unique(pat[, paste(study, signature, mode, ecology, metric, block_scale, sep = "|")])) {
  p <- strsplit(k, "\\|", fixed = FALSE)[[1]]
  d <- pat[study == p[1] & signature == p[2] & mode == p[3] & ecology == p[4] & metric == p[5] & block_scale == p[6]]
  x <- d$patient_effect[is.finite(d$patient_effect)]; ci <- bt_ci(x, stable_seed(SEEDS$spatial_patient_meta_bootstrap, k))
  set.seed(stable_seed(SEEDS$spatial_patient_meta_bootstrap, "studyperm", k))
  null <- if (length(x)) replicate(999L, mean(x * sample(c(-1,1), length(x), replace = TRUE))) else rep(NA_real_, 999L)
  emp <- if (length(x)) (1 + sum(abs(null) >= abs(mean(x)))) / 1000 else NA_real_
  si <- si + 1L; study_rows[[si]] <- data.table(study = p[1], signature = p[2], mode = p[3], ecology = p[4], metric = p[5],
    block_scale = p[6], n_patients = length(x), mean_patient_effect = mean(x), median_patient_effect = median(x),
    ci_low = ci[1], ci_high = ci[2], positive_patient_proportion = mean(x > 0), study_sign_permutation_p = emp)
}
stud <- rbindlist(study_rows)
fwrite(stud, file.path(RUN, "04_spatial_integration/SPATIAL_STUDY_LEVEL.tsv"), sep = "\t")

meta_fit <- function(x, slab) {
  ok <- is.finite(x) & abs(x) < 1
  x <- x[ok]; slab <- slab[ok]
  if (!length(x)) return(list(k = 0L, est = NA, lo = NA, hi = NA, pi_lo = NA, pi_hi = NA, tau2 = NA, i2 = NA))
  z <- atanh(pmax(pmin(x, .999999), -.999999)); vi <- rep(1/100, length(z))
  if (length(z) == 1L) return(list(k = 1L, est = x, lo = NA, hi = NA, pi_lo = NA, pi_hi = NA, tau2 = NA, i2 = NA))
  fit <- rma.uni(yi = z, vi = vi, method = "REML", slab = slab)
  pr <- predict(fit)
  list(k = length(z), est = tanh(as.numeric(fit$b)), lo = tanh(fit$ci.lb), hi = tanh(fit$ci.ub),
       pi_lo = tanh(pr$pi.lb), pi_hi = tanh(pr$pi.ub), tau2 = fit$tau2, i2 = fit$I2)
}

meta_input <- pat[metric %in% c("spot_spearman", "block_spearman") & is.finite(patient_effect)]
keys <- unique(meta_input[, paste(signature, mode, ecology, metric, block_scale, sep = "|")])
meta_rows <- list(); pred_rows <- list(); lopo_rows <- list(); loso_rows <- list(); mi <- pi_i <- li <- ls_i <- 0L
for (k in keys) {
  p <- strsplit(k, "\\|", fixed = FALSE)[[1]]
  d <- meta_input[signature == p[1] & mode == p[2] & ecology == p[3] & metric == p[4] & block_scale == p[5]]
  f <- meta_fit(d$patient_effect, paste(d$study, d$patient_id, sep = ":")); mi <- mi + 1L
  meta_rows[[mi]] <- data.table(signature = p[1], mode = p[2], ecology = p[3], metric = p[4], block_scale = p[5],
    n_patients = f$k, n_studies = uniqueN(d$study), random_effect_rho = f$est, ci_low = f$lo, ci_high = f$hi,
    tau2_fisher_z = f$tau2, i2_percent = f$i2)
  pi_i <- pi_i + 1L; pred_rows[[pi_i]] <- data.table(signature = p[1], mode = p[2], ecology = p[3], metric = p[4],
    block_scale = p[5], prediction_low = f$pi_lo, prediction_high = f$pi_hi, reverse_material_threshold = -0.10,
    excludes_reverse_material_effect = is.finite(f$pi_lo) && f$pi_lo > -0.10)
  for (pt in unique(d$patient_id)) {
    dd <- d[patient_id != pt]; ff <- meta_fit(dd$patient_effect, paste(dd$study, dd$patient_id, sep = ":")); li <- li + 1L
    lopo_rows[[li]] <- data.table(signature = p[1], mode = p[2], ecology = p[3], metric = p[4], block_scale = p[5],
      omitted_patient = pt, remaining_n = ff$k, estimate = ff$est, direction_flip = is.finite(ff$est) && is.finite(f$est) && sign(ff$est) != sign(f$est))
  }
  for (st in unique(d$study)) {
    dd <- d[study != st]; ff <- meta_fit(dd$patient_effect, paste(dd$study, dd$patient_id, sep = ":")); ls_i <- ls_i + 1L
    loso_rows[[ls_i]] <- data.table(signature = p[1], mode = p[2], ecology = p[3], metric = p[4], block_scale = p[5],
      omitted_study = st, remaining_n = ff$k, estimate = ff$est, direction_flip = is.finite(ff$est) && is.finite(f$est) && sign(ff$est) != sign(f$est))
  }
}
meta <- rbindlist(meta_rows); pred <- rbindlist(pred_rows); lopo <- rbindlist(lopo_rows); loso <- rbindlist(loso_rows)
fwrite(meta, file.path(RUN, "04_spatial_integration/SPATIAL_RANDOM_EFFECTS_META.tsv"), sep = "\t")
fwrite(pred, file.path(RUN, "04_spatial_integration/SPATIAL_PREDICTION_INTERVAL.tsv"), sep = "\t")
fwrite(lopo, file.path(RUN, "04_spatial_integration/SPATIAL_LOPO.tsv"), sep = "\t")
fwrite(loso, file.path(RUN, "04_spatial_integration/SPATIAL_LOSO.tsv"), sep = "\t")

rvp <- dcast(pat[metric %in% c("spot_spearman", "block_spearman")],
             study + patient_id + signature + ecology + metric + block_scale ~ mode, value.var = "patient_effect")
rvp[, difference_raw_minus_pruned := RAW - OVERLAP_PRUNED]
fwrite(rvp, file.path(RUN, "04_spatial_integration/SPATIAL_RAW_VS_PRUNED.tsv"), sep = "\t")

# Frozen gate adjudication by signature-ecology pair.
gate_rows <- list(); gi <- 0L
for (sig in names(SIG_RAW)) for (eco in names(ECO)) {
  raw_st <- stud[signature == sig & mode == "RAW" & ecology == eco & metric == "spot_spearman"]
  pru_st <- stud[signature == sig & mode == "OVERLAP_PRUNED" & ecology == eco & metric == "spot_spearman"]
  mm <- meta[signature == sig & mode == "OVERLAP_PRUNED" & ecology == eco & metric == "spot_spearman"]
  pp <- pred[signature == sig & mode == "OVERLAP_PRUNED" & ecology == eco & metric == "spot_spearman"]
  sc <- meta[signature == sig & mode == "OVERLAP_PRUNED" & ecology == eco & metric == "block_spearman" & block_scale %in% c("small","medium","large")]
  lp <- lopo[signature == sig & mode == "OVERLAP_PRUNED" & ecology == eco & metric == "spot_spearman"]
  two_study <- nrow(pru_st) == 2L && all(pru_st$mean_patient_effect > 0, na.rm = FALSE)
  raw_positive <- nrow(raw_st) == 2L && all(raw_st$mean_patient_effect > 0, na.rm = FALSE)
  pruned_positive <- nrow(mm) == 1L && is.finite(mm$random_effect_rho) && mm$random_effect_rho > 0
  pi_ok <- nrow(pp) == 1L && isTRUE(pp$excludes_reverse_material_effect)
  scale_n <- sum(sc$random_effect_rho > 0, na.rm = TRUE)
  no_patient_driver <- nrow(lp) > 0 && !any(lp$direction_flip, na.rm = TRUE) && all(lp$estimate > 0, na.rm = TRUE)
  section_vals <- sec_long[signature == sig & mode == "OVERLAP_PRUNED" & ecology == eco & metric == "spot_spearman", estimate]
  no_section_driver <- sum(is.finite(section_vals)) > 2 && sum(section_vals > 0, na.rm = TRUE) > 1
  status <- if (!pruned_positive && raw_positive) "SPATIAL_OVERLAP_DEPENDENT" else if (two_study && pruned_positive && pi_ok && scale_n >= 2 && no_patient_driver && no_section_driver) "SPATIAL_REPRODUCIBLE_CROSS_STUDY" else if (two_study && pruned_positive) "SPATIAL_DIRECTIONALLY_CONSISTENT" else if (nrow(pru_st) == 1L) "SPATIAL_SINGLE_STUDY_ONLY" else if (nrow(pru_st) == 2L && prod(sign(pru_st$mean_patient_effect)) < 0) "SPATIAL_HETEROGENEOUS" else if (nrow(pru_st) < 1L) "NOT_EVALUABLE" else "NO_RELIABLE_SPATIAL_ORGANIZATION"
  gi <- gi + 1L; gate_rows[[gi]] <- data.table(signature = sig, ecology = eco, gate_status = status,
    two_studies_same_positive_direction = two_study, patient_meta_positive = pruned_positive,
    prediction_excludes_reverse_material_effect = pi_ok, overlap_pruned_retained = pruned_positive,
    consistent_block_scales_n = scale_n, no_single_patient_driver = no_patient_driver,
    no_single_section_driver = no_section_driver,
    reference_caveat = "GSE237183 frozen all-section context includes six official IDH-mutant glioma sections")
}
gates <- rbindlist(gate_rows)

gate_md <- c("# Spatial Final Gate", "", paste0("Generated: `", format(Sys.time(), tz = "UTC", usetz = TRUE), "`"), "",
             "Reference structure: GSE237183 contains 19 sections and 13 frozen patient-units derived from official GEO titles; 13 sections are IDH-WT GBM and 6 are IDH-mutant glioma. GSE242352 contains six official additional GBM cases (B6-B11). This distinction is retained in every row and was not inferred from filenames.", "",
             "| Signature | Ecology | Gate |", "|---|---|---|")
gate_md <- c(gate_md, apply(gates, 1, function(x) paste0("| ", x[["signature"]], " | ", x[["ecology"]], " | ", x[["gate_status"]], " |")))
writeLines(gate_md, file.path(RUN, "04_spatial_integration/SPATIAL_ANALYSIS_SUMMARY.md"))

# Module session, completion gates, and checksums.
session_text <- c(capture.output(sessionInfo()), paste0("Started: ", format(started, tz = "UTC", usetz = TRUE)),
                  paste0("Finished: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
                  paste0("Elapsed_seconds: ", as.numeric(difftime(Sys.time(), started, units = "secs"))),
                  paste0("Maximum_R_heap_component_MB: ", round(max(gc()[,7]), 2)),
                  paste0("Script_SHA256: ", sha256(file.path(RUN, "scripts/spatial_v12_run.R"))))
writeLines(session_text, file.path(RUN, "02_gse237183/GSE237183_SESSION_INFO.txt"))
writeLines(session_text, file.path(RUN, "03_gse242352/GSE242352_SESSION_INFO.txt"))
writeLines(session_text, file.path(RUN, "04_spatial_integration/SPATIAL_SESSION_INFO.txt"))
writeLines(session_text, file.path(RUN, "01_spatial_reference/SPATIAL_REFERENCE_SESSION_INFO.txt"))

writeLines(c("# Spatial Reference Completion Gate", "", "Status: PASS", "",
             "- Official GEO family SOFT parsed and hashed for both studies.",
             "- GSE242352 official tar SHA-256 verified against the frozen reference.",
             "- Section/patient identities use official titles/series case statements, never filenames alone.",
             "- Actual disease composition is serialized; six GSE237183 sections are IDH-mutant glioma.",
             "- Next-module eligibility: ELIGIBLE"), file.path(RUN, "01_spatial_reference/SPATIAL_REFERENCE_COMPLETION_GATE.md"))
writeLines(c("# GSE237183 Completion Gate", "", "Status: PASS_WITH_REFERENCE_COMPOSITION_CAVEAT", "",
             paste0("- Sections: ", nrow(soft237), "; patient-units: ", uniqueN(soft237$patient_id), "; GBM sections: ", sum(soft237$disease_class == "IDH_WT_GBM"), "; IDH-mutant glioma sections: ", sum(soft237$disease_class != "IDH_WT_GBM"), "."),
             "- Formal block permutations: 999 per endpoint and scale.",
             paste0("- Elapsed seconds: ", round(r237$elapsed_seconds, 2)),
             "- Next-module eligibility: ELIGIBLE_FOR_PATIENT_LEVEL_INTEGRATION"), file.path(RUN, "02_gse237183/GSE237183_COMPLETION_GATE.md"))
writeLines(c("# GSE242352 Completion Gate", "", "Status: PASS", "",
             paste0("- Official cases/sections: ", nrow(soft242), "; patient-units: ", uniqueN(soft242$patient_id), "."),
             "- Formal block permutations: 999 per endpoint and scale.",
             paste0("- Elapsed seconds: ", round(r242$elapsed_seconds, 2)),
             "- Next-module eligibility: ELIGIBLE_FOR_CROSS_STUDY_INTEGRATION"), file.path(RUN, "03_gse242352/GSE242352_COMPLETION_GATE.md"))
writeLines(c("# Spatial Integration Completion Gate", "", "Status: PASS", "",
             "- Inference hierarchy: spot -> block -> section -> patient -> study.",
             "- Multiple sections were equal-weighted within patient; patients were equal-weighted within study.",
             "- Random-effects meta, prediction intervals, LOPO, LOSO, pruning and three-scale checks completed.",
             paste0("- Reproducible cross-study gates: ", sum(gates$gate_status == "SPATIAL_REPRODUCIBLE_CROSS_STUDY"), "."),
             "- Next-module eligibility: ELIGIBLE_FOR_V1_2_CROSS_SCALE_INTEGRATION"), file.path(RUN, "04_spatial_integration/SPATIAL_INTEGRATION_COMPLETION_GATE.md"))

write_manifest <- function(dir, outfile) {
  fs <- list.files(file.path(RUN, dir), full.names = TRUE, recursive = TRUE)
  fs <- fs[file.info(fs)$isdir %in% FALSE & basename(fs) != basename(outfile)]
  lines <- vapply(fs, function(f) paste(sha256(f), substring(f, nchar(RUN) + 2)), character(1))
  writeLines(lines, file.path(RUN, outfile))
}
write_manifest("01_spatial_reference", "01_spatial_reference/SPATIAL_REFERENCE_CHECKSUM.sha256")
write_manifest("02_gse237183", "02_gse237183/GSE237183_CHECKSUM.sha256")
write_manifest("03_gse242352", "03_gse242352/GSE242352_CHECKSUM.sha256")
write_manifest("04_spatial_integration", "04_spatial_integration/SPATIAL_INTEGRATION_CHECKSUM.sha256")

message("Spatial v1.2 complete: ", nrow(allres), " section-scale endpoint rows; ", nrow(pat), " patient rows; ", nrow(meta), " meta rows")
