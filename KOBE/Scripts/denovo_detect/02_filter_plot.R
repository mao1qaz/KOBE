#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11L) {
  stop(paste(
    "Usage: Rscript 02_filter_plot.R INPUT PREFIX REGION SAMPLE_TABLE SPECIES",
    "DNA_MIN DNA_FDR_CUTOFF RNA_MIN RNA_FDR_CUTOFF TOTAL_ERROR_RATE MIN_SUPPORT"
  ))
}

infile <- args[1]
prefix <- args[2]
region <- args[3]
sample_table <- args[4]
species_id <- args[5]
dna_min <- as.integer(args[6])
dna_fdr_cutoff <- as.numeric(args[7])
rna_min <- as.integer(args[8])
rna_fdr_cutoff <- as.numeric(args[9])
total_error_rate <- as.numeric(args[10])
min_support <- as.integer(args[11])

if (!(region %in% c("genic", "intergenic"))) stop("Invalid region: ", region)
if (is.na(dna_min) || dna_min < 1L || is.na(rna_min) || rna_min < 1L) stop("Coverage thresholds must be positive integers")
if (!is.finite(dna_fdr_cutoff) || dna_fdr_cutoff <= 0 || dna_fdr_cutoff >= 1) stop("Invalid DNA FDR cutoff")
if (!is.finite(rna_fdr_cutoff) || rna_fdr_cutoff <= 0 || rna_fdr_cutoff >= 1) stop("Invalid RNA FDR cutoff")
if (!is.finite(total_error_rate) || total_error_rate <= 0 || total_error_rate >= 1) stop("Invalid total sequencing error rate")
if (is.na(min_support) || min_support < 1L) stop("min_support must be a positive integer")
if (!file.exists(sample_table)) stop("Sample table not found: ", sample_table)

# Read and validate the sample manifest. DNA/RNA matching is performed by sample_name.
meta <- fread(sample_table, sep = "\t", header = TRUE, data.table = TRUE, check.names = FALSE, showProgress = FALSE)
required_meta <- c("sample_id", "sample_name", "species", "type")
missing_meta <- setdiff(required_meta, names(meta))
if (length(missing_meta)) stop("Sample table is missing columns: ", paste(missing_meta, collapse = ", "))
meta[, sample_id := trimws(as.character(sample_id))]
meta[, sample_name := trimws(as.character(sample_name))]
meta[, species := trimws(as.character(species))]
meta[, type := toupper(trimws(as.character(type)))]
meta <- meta[species == species_id]
if (!nrow(meta)) stop("No samples found for species: ", species_id)
if (anyNA(meta$sample_id) || any(meta$sample_id == "") || anyNA(meta$sample_name) || any(meta$sample_name == "")) {
  stop("sample_id and sample_name must be non-empty for every sample")
}
if (anyDuplicated(meta$sample_id)) stop("sample_id must be unique within the species: ", species_id)
if (any(!meta$type %in% c("DNA", "RNA"))) stop("Sample type must be DNA or RNA for species: ", species_id)

dna_libs <- meta[type == "DNA", sample_id]
rna_libs <- meta[type == "RNA", sample_id]
if (!length(dna_libs) || !length(rna_libs)) stop("At least one DNA and one RNA library are required for species: ", species_id)

rna_groups <- unique(meta[type == "RNA", sample_name])
for (sn in rna_groups) {
  group_meta <- meta[sample_name == sn]
  if (!any(group_meta$type == "DNA")) {
    stop("RNA sample_name has no matched DNA library: ", sn)
  }
}

# Read the annotated count table. Each library retains its own reference/ALT counts.
df <- fread(infile, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)
required_base <- c("chr", "pos", "ref", "alt", "anno", "gene_id", "transcript", "cds_pos", "mut", "AA_pos", "AA_change", "region")
required_counts <- as.vector(rbind(paste0("R_", c(dna_libs, rna_libs)), paste0("A_", c(dna_libs, rna_libs))))
missing_cols <- setdiff(c(required_base, required_counts), names(df))
if (length(missing_cols)) stop("Input count table is missing columns: ", paste(missing_cols, collapse = ", "))
if (nrow(df) && any(as.character(df$region) != region)) stop("Input contains records from another region")

# The legacy binomial error model is retained. ALT reads are tested as X >= ALT.
eps0 <- (total_error_rate / 3) / (1 - 2 * total_error_rate / 3)
calc_p <- function(alt, cov) {
  p <- pbinom(alt - 1, size = cov, prob = eps0, lower.tail = FALSE)
  p[!is.finite(p)] <- 1
  p
}

all_libs <- c(dna_libs, rna_libs)
for (sid in all_libs) {
  r_col <- paste0("R_", sid)
  a_col <- paste0("A_", sid)
  c_col <- paste0("C_", sid)
  l_col <- paste0("L_", sid)
  p_col <- paste0("perr_", sid)
  q_col <- paste0("padj_", sid)
  df[[r_col]] <- as.numeric(df[[r_col]])
  df[[a_col]] <- as.numeric(df[[a_col]])
  if (anyNA(df[[r_col]]) || anyNA(df[[a_col]])) stop("Missing or non-numeric read counts for sample_id: ", sid)
  if (any(df[[r_col]] < 0 | df[[a_col]] < 0)) stop("Negative read counts for sample_id: ", sid)
  df[[c_col]] <- df[[r_col]] + df[[a_col]]
  df[[l_col]] <- fifelse(df[[c_col]] > 0, df[[a_col]] / df[[c_col]], NA_real_)
  df[[p_col]] <- NA_real_
  df[[q_col]] <- NA_real_
}

# Each sample_name is an independent DNA-background group.
# For each group, ALL its DNA libraries must have coverage >= dna_min.
# Each DNA library is tested separately and BH-adjusted over that group's DNA-coverage universe.
group_pass_mat <- matrix(FALSE, nrow = nrow(df), ncol = length(rna_groups), dimnames = list(NULL, rna_groups))
rna_pass_mat <- matrix(FALSE, nrow = nrow(df), ncol = length(rna_libs), dimnames = list(NULL, rna_libs))
dna_summary_rows <- list()
rna_summary_rows <- list()

for (g_idx in seq_along(rna_groups)) {
  sn <- rna_groups[g_idx]
  group_meta <- meta[sample_name == sn]
  group_dna <- group_meta[type == "DNA", sample_id]
  group_rna <- group_meta[type == "RNA", sample_id]

  dna_cov_mat <- do.call(cbind, lapply(group_dna, function(sid) df[[paste0("C_", sid)]]))
  if (is.null(dim(dna_cov_mat))) dna_cov_mat <- matrix(dna_cov_mat, ncol = 1L)
  dna_cov_ok <- if (nrow(df)) rowSums(dna_cov_mat >= dna_min) == length(group_dna) else logical(0)
  dna_idx <- which(dna_cov_ok)
  dna_lib_pass_mat <- matrix(FALSE, nrow = length(dna_idx), ncol = length(group_dna), dimnames = list(NULL, group_dna))
  dna_lib_pass_counts <- integer(length(group_dna))

  if (length(dna_idx)) {
    for (d_idx in seq_along(group_dna)) {
      sid <- group_dna[d_idx]
      cov <- df[[paste0("C_", sid)]][dna_idx]
      alt <- df[[paste0("A_", sid)]][dna_idx]
      p_values <- calc_p(alt, cov)
      q_values <- p.adjust(p_values, method = "BH")
      df[[paste0("perr_", sid)]][dna_idx] <- p_values
      df[[paste0("padj_", sid)]][dna_idx] <- q_values
      dna_lib_pass_mat[, d_idx] <- q_values > dna_fdr_cutoff
      dna_lib_pass_counts[d_idx] <- sum(q_values > dna_fdr_cutoff)
    }
    group_dna_pass <- logical(nrow(df))
    group_dna_pass[dna_idx] <- rowSums(dna_lib_pass_mat) == length(group_dna)
  } else {
    group_dna_pass <- logical(nrow(df))
  }
  group_pass_mat[, g_idx] <- group_dna_pass
  n_group_dna_pass <- sum(group_dna_pass)

  for (d_idx in seq_along(group_dna)) {
    dna_summary_rows[[length(dna_summary_rows) + 1L]] <- data.table(
      species = species_id, region = region, type = "DNA", sample_id = group_dna[d_idx],
      sample_name = sn, matched_dna_libraries = paste(group_dna, collapse = ";"),
      dna_coverage_test_sites = length(dna_idx),
      dna_library_fdr_pass_sites = dna_lib_pass_counts[d_idx],
      dna_group_pass_sites = n_group_dna_pass,
      rna_test_sites = NA_integer_, rna_fdr_pass_sites = NA_integer_,
      final_supported_sites_for_library = NA_integer_
    )
  }

  # Only after the complete DNA group passes do we test each matched RNA library.
  for (sid in group_rna) {
    rna_cov_ok <- df[[paste0("C_", sid)]] >= rna_min
    eligible <- group_dna_pass & rna_cov_ok
    rna_idx <- which(eligible)
    rna_pass <- logical(nrow(df))
    if (length(rna_idx)) {
      cov <- df[[paste0("C_", sid)]][rna_idx]
      alt <- df[[paste0("A_", sid)]][rna_idx]
      p_values <- calc_p(alt, cov)
      q_values <- p.adjust(p_values, method = "BH")
      df[[paste0("perr_", sid)]][rna_idx] <- p_values
      df[[paste0("padj_", sid)]][rna_idx] <- q_values
      rna_pass[rna_idx] <- q_values < rna_fdr_cutoff
    }
    rna_col_idx <- match(sid, rna_libs)
    # Libraries sharing sample_name remain separate, independent columns.
    rna_pass_mat[, rna_col_idx] <- rna_pass
    rna_summary_rows[[length(rna_summary_rows) + 1L]] <- data.table(
      species = species_id, region = region, type = "RNA", sample_id = sid,
      sample_name = sn, matched_dna_libraries = paste(group_dna, collapse = ";"),
      dna_coverage_test_sites = length(dna_idx),
      dna_library_fdr_pass_sites = NA_integer_,
      dna_group_pass_sites = n_group_dna_pass,
      rna_test_sites = length(rna_idx), rna_fdr_pass_sites = sum(rna_pass),
      final_supported_sites_for_library = NA_integer_
    )
  }
}

# RNA pass matrix is filled per sample_id during the matched DNA-then-RNA loop.

# Summarize DNA groups that pass for each site and independent RNA library support.
if (nrow(df)) {
  df[, DNA_pass_sample_names := vapply(seq_len(.N), function(i) paste(rna_groups[group_pass_mat[i, ]], collapse = ";"), character(1))]
  df[, RNA_pass_samples := vapply(seq_len(.N), function(i) paste(rna_libs[rna_pass_mat[i, ]], collapse = ";"), character(1))]
  df[, n_RNA_pass := rowSums(rna_pass_mat)]
  df[, edit_level_pass_samples := vapply(seq_len(.N), function(i) {
    selected <- rna_libs[rna_pass_mat[i, ]]
    if (!length(selected)) return("")
    values <- vapply(selected, function(sid) as.numeric(df[[paste0("L_", sid)]][i]), numeric(1))
    paste(paste0(selected, ":", formatC(values, format = "f", digits = 4)), collapse = ";")
  }, character(1))]
  if (region == "genic") {
  df[, anno_group := as.character(anno)]
} else {
  df[, anno_group := "Intergenic"]
}
} else {
  df[, `:=`(DNA_pass_sample_names = character(), RNA_pass_samples = character(), n_RNA_pass = integer(), edit_level_pass_samples = character(), anno_group = character())]
}

# Any RNA library passing is reported separately from the stricter final support requirement.
df_any_rna <- df[n_RNA_pass >= 1L]
df_supported <- df[n_RNA_pass >= min_support]
selected_mut <- if (region == "genic") "A>G" else c("A>G", "T>C")
df_selected <- df_supported[mut %in% selected_mut]

# Update per-RNA-library final-support counts in the sample-level audit table.
for (i in seq_along(rna_libs)) {
  sid <- rna_libs[i]
  supported_for_lib <- if (nrow(df)) sum(rna_pass_mat[, i] & (rowSums(rna_pass_mat) >= min_support)) else 0L
  for (j in seq_along(rna_summary_rows)) {
    if (rna_summary_rows[[j]]$sample_id[1] == sid) {
      rna_summary_rows[[j]][, final_supported_sites_for_library := supported_for_lib]
      break
    }
  }
}

# Write detailed per-site statistics, including per-library raw P and BH-FDR columns.
fwrite(df, paste0(prefix, ".with_stat.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(df_any_rna, paste0(prefix, ".fdr05.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(df_supported, paste0(prefix, ".supported.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(df_selected, paste0(prefix, ".selected.tsv"), sep = "\t", quote = FALSE, na = "NA")

sample_summary <- rbindlist(c(dna_summary_rows, rna_summary_rows), fill = TRUE)
setorder(sample_summary, type, sample_name, sample_id)
fwrite(sample_summary, paste0(prefix, ".sample_test_summary.tsv"), sep = "\t", quote = FALSE, na = "NA")

# Region-level summary counts groups and independent library tests separately.
dna_group_site_pairs <- if (ncol(group_pass_mat)) sum(group_pass_mat) else 0L
sites_with_any_dna_group_pass <- if (nrow(group_pass_mat)) sum(rowSums(group_pass_mat) > 0L) else 0L
rna_test_site_library_pairs <- sum(vapply(rna_summary_rows, function(x) as.numeric(x$rna_test_sites[1]), numeric(1)))
rna_fdr_pass_site_library_pairs <- sum(vapply(rna_summary_rows, function(x) as.numeric(x$rna_fdr_pass_sites[1]), numeric(1)))
summary_dt <- data.table(
  species = species_id, region = region, candidate_sites = nrow(df),
  sample_name_groups = length(rna_groups),
  dna_group_pass_site_pairs = dna_group_site_pairs,
  sites_with_any_dna_group_pass = sites_with_any_dna_group_pass,
  rna_test_site_library_pairs = rna_test_site_library_pairs,
  rna_fdr_pass_site_library_pairs = rna_fdr_pass_site_library_pairs,
  sites_with_any_rna_pass = nrow(df_any_rna),
  sites_meeting_min_rna_support = nrow(df_supported),
  final_selected_sites = nrow(df_selected),
  dna_min_coverage = dna_min, dna_fdr_must_be_greater_than = dna_fdr_cutoff,
  rna_min_coverage = rna_min, rna_fdr_must_be_less_than = rna_fdr_cutoff,
  total_error_rate = total_error_rate, min_supporting_rna_sample_ids = min_support
)
fwrite(summary_dt, paste0(prefix, ".filter_summary.tsv"), sep = "\t", quote = FALSE, na = "NA")

# Plot mutation classes among sites meeting the final independent-RNA support requirement.
mut_order <- c("A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G")
mut_count <- as.integer(table(factor(df_supported$mut, levels = mut_order)))
bar_cols <- rep("aquamarine2", length(mut_order))
bar_cols[mut_order %in% c("A>G", "T>C")] <- "#cd2626"
pdf(paste0(prefix, ".mutation_types.pdf"), width = 8, height = 5.5, useDingbats = FALSE)
par(mar = c(6, 5, 3, 1), cex.axis = 0.9, cex.lab = 1.05)
ymax <- max(1, mut_count) * 1.18
bp <- barplot(mut_count, col = bar_cols, border = "black", names.arg = mut_order, las = 2,
              ylab = "Number of sites", main = paste("Mutation types among supported", region, "sites"), ylim = c(0, ymax))
text(bp, mut_count + ymax * 0.035, labels = mut_count, cex = 0.75, xpd = NA)
dev.off()

cat("Species:", species_id, "\nRegion:", region,
    "\nCandidate sites:", nrow(df),
    "\nSites with any passing DNA sample_name group:", sites_with_any_dna_group_pass,
    "\nSites with any RNA library passing:", nrow(df_any_rna),
    "\nSites meeting minimum RNA sample_id support:", nrow(df_supported),
    "\nFinal selected sites:", nrow(df_selected),
    "\nDNA coverage >=", dna_min, "; DNA BH-FDR >", dna_fdr_cutoff,
    "\nRNA coverage >=", rna_min, "; RNA BH-FDR <", rna_fdr_cutoff,
    "\nMinimum independent RNA sample_id support:", min_support,
    "\nTotal sequencing error rate:", total_error_rate, "\n")
