#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

# Resolve helper scripts relative to this R script.
script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_file_arg) != 1L) {
    stop("Cannot determine the path of run_pairwise_tests.R.")
}
script_file <- normalizePath(
    sub("^--file=", "", script_file_arg[1]),
    winslash = "/",
    mustWork = TRUE
)
SCRIPT_DIR <- dirname(script_file)

# Statistical parameters
DEFAULT_MIN_COVERAGE <- 3
TTTEST_MIN <- 2
FISHER_MIN <- 1
REDIT_MIN <- 1
DEFAULT_FDR_CUTOFF <- 0.05
FDR_CUTOFF <- DEFAULT_FDR_CUTOFF
TOP_LABELS <- 20

# Usage
print_usage <- function() {
    cat("\nRNA editing pairwise hypothesis testing\n\n",
        "Required arguments:\n",
        "  --input FILE              Input editing count TSV\n",
        "  --group STRING            Define a group (repeatable)\n",
        "  --comparisons STRING     Specify pairwise comparisons\n",
        "  --outdir DIR              Output parent directory\n\n",
        "Optional arguments:\n",
        "  --min_coverage INT        Minimum A+G coverage (default: 3)\n",
        "  --fdr_cutoff FLOAT        BH-FDR significance cutoff (default: 0.05)\n\n",
        "Example:\n\n",
        "Rscript /path/to/run_pairwise_tests.R \\\n",
        "    --input /path/to/editing_counts.tsv \\\n",
        "    --group N=N1,N2,N3 \\\n",
        "    --group F=F1,F2,F3 \\\n",
        "    --comparisons N:F \\\n",
        "    --min_coverage 3 \\\n",
        "    --fdr_cutoff 0.05 \\\n",
        "    --outdir /path/to/diff_result\n\n", sep = "")
}

# Parse command-line arguments
parse_args <- function(args) {
    opt <- list(input = NULL, outdir = NULL, min_coverage = DEFAULT_MIN_COVERAGE,
                fdr_cutoff = DEFAULT_FDR_CUTOFF, groups = list(), comparisons = NULL)
    i <- 1L
    while (i <= length(args)) {
        key <- args[i]
        if (key == "--help") { print_usage(); quit(status = 0) }
        if (i == length(args)) stop("Missing value for argument: ", key)
        value <- args[i + 1L]
        if (key == "--input") {
            opt$input <- value
        } else if (key == "--outdir") {
            opt$outdir <- value
        } else if (key == "--min_coverage") {
            opt$min_coverage <- suppressWarnings(as.numeric(value))
        } else if (key == "--fdr_cutoff") {
            opt$fdr_cutoff <- suppressWarnings(as.numeric(value))
        } else if (key == "--group") {
            eq <- regexpr("=", value, fixed = TRUE)[1]
            if (eq < 2 || eq == nchar(value)) stop("Invalid group format: ", value)
            group_name <- substr(value, 1, eq - 1)
            sample_string <- substr(value, eq + 1, nchar(value))
            if (grepl("[,:]", group_name)) stop("Group names cannot contain ':' or ','.")
            if (group_name %in% names(opt$groups)) stop("Duplicate group name: ", group_name)
            group_samples <- trimws(strsplit(sample_string, ",", fixed = TRUE)[[1]])
            if (any(group_samples == "") || anyDuplicated(group_samples))
                stop("Invalid or duplicated samples in group: ", group_name)
            opt$groups[[group_name]] <- group_samples
        } else if (key == "--comparisons") {
            if (!is.null(opt$comparisons)) stop("Specify --comparisons only once.")
            opt$comparisons <- value
        } else {
            stop("Unknown argument: ", key)
        }
        i <- i + 2L
    }

    if (is.null(opt$input) || opt$input == "") {
        cat("\nERROR: --input is required.\n"); print_usage(); quit(status = 1)
    }
    if (is.null(opt$outdir) || !nzchar(opt$outdir))
        stop("--outdir is required.")
    if (!file.exists(opt$input)) stop("Input file not found: ", opt$input)
    if (length(opt$groups) < 2) stop("At least two --group arguments are required.")
    if (is.null(opt$comparisons)) stop("--comparisons is required.")
    if (length(opt$min_coverage) != 1 || !is.finite(opt$min_coverage) ||
        opt$min_coverage < 1 || opt$min_coverage != floor(opt$min_coverage))
        stop("--min_coverage must be a positive integer.")
    if (length(opt$fdr_cutoff) != 1 || !is.finite(opt$fdr_cutoff) ||
        opt$fdr_cutoff <= 0 || opt$fdr_cutoff > 1)
        stop("--fdr_cutoff must be greater than 0 and less than or equal to 1.")
    all_samples <- unlist(opt$groups, use.names = FALSE)
    if (anyDuplicated(all_samples)) stop("A sample cannot belong to multiple groups.")
    pairs <- strsplit(opt$comparisons, ",", fixed = TRUE)[[1]]
    comparisons <- lapply(pairs, function(x) {
        z <- strsplit(trimws(x), ":", fixed = TRUE)[[1]]
        if (length(z) != 2 || any(z == "")) stop("Invalid comparison: ", x)
        if (z[1] == z[2]) stop("A group cannot be compared with itself: ", x)
        if (!all(z %in% names(opt$groups))) stop("Comparison contains an undefined group: ", x)
        z
    })
    comparison_names <- vapply(comparisons, function(x) paste(x, collapse = "_vs_"), character(1))
    if (anyDuplicated(comparison_names)) stop("Duplicate comparisons were specified.")
    opt$comparisons <- comparisons
    opt
}

# General helpers
get_named_value <- function(x, name) {
    if (is.null(x)) return(NULL)
    if (!is.null(names(x)) && name %in% names(x)) return(x[[name]])
    if (is.data.frame(x) && name %in% names(x)) return(x[[name]][1])
    if (is.matrix(x) && !is.null(rownames(x)) && name %in% rownames(x)) return(x[name, 1])
    NULL
}

get_numeric <- function(x, name = NULL) {
    if (!is.null(name)) x <- get_named_value(x, name)
    if (is.null(x) || length(x) == 0) return(NA_real_)
    value <- suppressWarnings(as.numeric(x[1]))
    if (length(value) == 0 || !is.finite(value)) return(NA_real_)
    value
}

get_logical <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x[1])) return(NA)
    as.logical(x[1])
}

bh_adjust <- function(p) {
    result <- rep(NA_real_, length(p)); keep <- !is.na(p) & is.finite(p)
    if (any(keep)) result[keep] <- p.adjust(p[keep], method = "BH")
    result
}

run_test <- function(fun) {
    tryCatch(list(p = fun(), error = NA_character_),
             error = function(e) list(p = NA_real_, error = conditionMessage(e)))
}

test_status <- function(p, err, eligible) {
    if (!eligible) return("insufficient_valid_samples")
    if (!is.na(err) && nzchar(err)) return("error")
    if (is.na(p) || !is.finite(p)) return("p_value_unavailable")
    "tested"
}

# Load REDIT functions
load_redit <- function() {
    llr_file <- file.path(SCRIPT_DIR, "REDIT_LLR.R")
    reg_file <- file.path(SCRIPT_DIR, "REDIT_regression.R")
    if (!file.exists(llr_file)) stop("Missing REDIT-LLR script: ", llr_file)
    if (!file.exists(reg_file)) stop("Missing REDIT-Regression script: ", reg_file)
    LLR_env <- new.env(parent = globalenv()); REG_env <- new.env(parent = globalenv())
    sys.source(llr_file, envir = LLR_env); sys.source(reg_file, envir = REG_env)
    if (!is.function(LLR_env$REDIT_LLR)) stop("REDIT_LLR() not found.")
    if (!is.function(REG_env$REDIT_regression)) stop("REDIT_regression() not found.")
    list(LLR = LLR_env$REDIT_LLR, Regression = REG_env$REDIT_regression)
}

# Read and validate input counts
read_counts <- function(input_file, groups) {
    dat <- fread(input_file, sep = "\t", check.names = FALSE, data.table = FALSE)
    if (!"ID" %in% names(dat)) stop("Missing ID column.")
    if (anyNA(dat$ID) || any(dat$ID == "")) stop("Missing site IDs.")
    if (anyDuplicated(dat$ID)) stop("Duplicate site IDs.")
    a_cols <- grep("_A_reads$", names(dat), value = TRUE)
    g_cols <- grep("_G_reads$", names(dat), value = TRUE)
    a_samples <- sub("_A_reads$", "", a_cols); g_samples <- sub("_G_reads$", "", g_cols)
    if (length(a_samples) == 0 || length(g_samples) == 0) stop("No A_reads/G_reads columns found.")
    if (!setequal(a_samples, g_samples)) stop("A_reads and G_reads sample columns do not match.")
    if (anyDuplicated(a_samples) || anyDuplicated(g_samples)) stop("Duplicate sample column names.")
    requested <- unique(unlist(groups, use.names = FALSE))
    missing_samples <- setdiff(requested, a_samples)
    if (length(missing_samples) > 0)
        stop("Samples not found in input matrix: ", paste(missing_samples, collapse = ", "))
    for (sample in requested) for (suffix in c("_A_reads", "_G_reads")) {
        col <- paste0(sample, suffix); dat[[col]] <- suppressWarnings(as.numeric(dat[[col]]))
        x <- dat[[col]]
        if (anyNA(x) || any(!is.finite(x)) || any(x < 0) || any(x != floor(x)))
            stop("Invalid read counts in column: ", col)
    }
    cat("Input:", input_file, "\nTotal sites:", nrow(dat),
        "\nDetected RNA samples:", length(a_samples), "\n")
    dat
}

# Plot one statistical-method panel
plot_panel <- function(df, fdr_col, method_name, g1, g2, x_col, y_col, level_label) {
    valid <- is.finite(df[[x_col]]) & is.finite(df[[y_col]])
    plot_df <- df[valid, , drop = FALSE]
    xlab <- paste0(g1, " ", level_label); ylab <- paste0(g2, " ", level_label)
    if (nrow(plot_df) == 0) {
        plot(NA, xlim = c(0, 1), ylim = c(0, 1), asp = 1, xlab = xlab, ylab = ylab,
             main = method_name, axes = FALSE)
        axis(1, at = seq(0, 1, 0.2)); axis(2, at = seq(0, 1, 0.2), las = 1)
        box(bty = "l"); text(0.5, 0.5, "No sites with valid editing levels")
        return(invisible(NULL))
    }
    fdr <- plot_df[[fdr_col]]; significant <- !is.na(fdr) & fdr < FDR_CUTOFF
    selected <- plot_df[significant, , drop = FALSE]; n_sig <- nrow(selected)
    plot(plot_df[[x_col]], plot_df[[y_col]], xlim = c(0, 1), ylim = c(0, 1), asp = 1,
         xlab = xlab, ylab = ylab, pch = 16, col = "grey75", cex = 0.7, axes = FALSE,
         main = paste0(method_name, "\nFDR < ", FDR_CUTOFF, ", n = ", n_sig))
    axis(1, at = seq(0, 1, 0.2)); axis(2, at = seq(0, 1, 0.2), las = 1)
    abline(0, 1, lty = 2, col = "grey50"); box(bty = "l")
    if (n_sig > 0) {
        points(selected[[x_col]], selected[[y_col]], pch = 16, col = "red", cex = 1)
        selected <- selected[order(selected[[fdr_col]], selected$ID), , drop = FALSE]
        labels_df <- head(selected, TOP_LABELS)
        text(labels_df[[x_col]], labels_df[[y_col]], labels = labels_df$ID,
             pos = 3, offset = 0.35, cex = 0.45, col = "black", xpd = NA)
    }
    legend("topleft", legend = c("All sites", paste0("FDR < ", FDR_CUTOFF)), col = c("grey75", "red"),
           pch = 16, bty = "n", cex = 0.65)
    invisible(NULL)
}

# Generate a four-panel PDF: mean or pooled editing levels
plot_pair <- function(df, g1, g2, result_dir, level_type = c("mean", "pooled")) {
    level_type <- match.arg(level_type); comparison_name <- paste0(g1, "_vs_", g2)
    if (level_type == "mean") {
        x_col <- "mean_group1"; y_col <- "mean_group2"
        level_label <- "mean editing level"; plot_title <- "Mean editing levels"
        file_suffix <- "_scatter.pdf"
    } else {
        x_col <- "pooled_group1"; y_col <- "pooled_group2"
        level_label <- "pooled editing level"; plot_title <- "Pooled editing levels"
        file_suffix <- "_pooled_scatter.pdf"
    }
    pdf_file <- file.path(result_dir, paste0(comparison_name, file_suffix))
    pdf(pdf_file, width = 12, height = 9, onefile = TRUE)
    on.exit(dev.off(), add = TRUE)
    par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.5, 1.5), oma = c(0, 0, 1.5, 0),
        bty = "l", pty = "s")
    plot_panel(df, "ttest_FDR", "Student's t-test", g1, g2, x_col, y_col, level_label)
    plot_panel(df, "fisher_FDR", "Fisher's exact test", g1, g2, x_col, y_col, level_label)
    plot_panel(df, "LLR_FDR", "REDIT-LLR", g1, g2, x_col, y_col, level_label)
    plot_panel(df, "Regression_FDR", "REDIT-Regression", g1, g2, x_col, y_col, level_label)
    mtext(paste0(comparison_name, " | ", plot_title), outer = TRUE, side = 3,
          line = 0.3, font = 2, cex = 1.1)
    cat("Scatter plot saved:", pdf_file, "\n")
    invisible(pdf_file)
}

# Run one pairwise comparison
run_pair <- function(counts, g1, g2, samples1, samples2, min_coverage, redit, outdir) {
    cat("\n========================================\nComparison:", g1, "vs", g2,
        "\n========================================\n")
    cat("Group 1:", paste(samples1, collapse = ", "), "\n")
    cat("Group 2:", paste(samples2, collapse = ", "), "\n")
    cat("Minimum coverage:", min_coverage, "\n")
    n_sites <- nrow(counts); results <- vector("list", n_sites)

    for (i in seq_len(n_sites)) {
        site <- counts$ID[i]
        a1 <- vapply(samples1, function(s) counts[[paste0(s, "_A_reads")]][i], numeric(1))
        g1_reads <- vapply(samples1, function(s) counts[[paste0(s, "_G_reads")]][i], numeric(1))
        a2 <- vapply(samples2, function(s) counts[[paste0(s, "_A_reads")]][i], numeric(1))
        g2_reads <- vapply(samples2, function(s) counts[[paste0(s, "_G_reads")]][i], numeric(1))
        cov1 <- a1 + g1_reads; cov2 <- a2 + g2_reads
        keep1 <- cov1 >= min_coverage; keep2 <- cov2 >= min_coverage
        a1 <- a1[keep1]; g1_reads <- g1_reads[keep1]; cov1 <- cov1[keep1]
        a2 <- a2[keep2]; g2_reads <- g2_reads[keep2]; cov2 <- cov2[keep2]
        n1 <- length(cov1); n2 <- length(cov2)
        valid_samples1 <- samples1[keep1]; valid_samples2 <- samples2[keep2]
        mean1 <- if (n1 > 0) mean(g1_reads / cov1) else NA_real_
        mean2 <- if (n2 > 0) mean(g2_reads / cov2) else NA_real_
        pooled1 <- if (n1 > 0) sum(g1_reads) / sum(cov1) else NA_real_
        pooled2 <- if (n2 > 0) sum(g2_reads) / sum(cov2) else NA_real_

        ttest_p <- fisher_p <- LLR_p <- Regression_p <- NA_real_
        Regression_effect <- Regression_intercept <- Regression_sigma <- NA_real_
        Regression_converged <- NA
        ttest_error <- fisher_error <- LLR_error <- Regression_error <- NA_character_

        # Student's t-test (equal variance)
        ttest_eligible <- n1 >= TTTEST_MIN && n2 >= TTTEST_MIN
        if (ttest_eligible) {
            ans <- run_test(function() t.test(g1_reads / cov1, g2_reads / cov2,
                                               alternative = "two.sided", var.equal = TRUE)$p.value)
            ttest_p <- ans$p; ttest_error <- ans$error
        }

        # Fisher's exact test
        fisher_eligible <- n1 >= FISHER_MIN && n2 >= FISHER_MIN
        if (fisher_eligible) {
            ans <- run_test(function() {
                tab <- matrix(c(sum(g1_reads), sum(a1), sum(g2_reads), sum(a2)), nrow = 2, byrow = TRUE)
                fisher.test(tab)$p.value
            })
            fisher_p <- ans$p; fisher_error <- ans$error
        }

        # REDIT-LLR and REDIT-Regression
        redit_eligible <- n1 >= REDIT_MIN && n2 >= REDIT_MIN
        if (redit_eligible) {
            redit_groups <- c(rep(g1, n1), rep(g2, n2))
            redit_data <- rbind(edited = c(g1_reads, g2_reads), unedited = c(a1, a2))
            llr_ans <- run_test(function() {
                llr <- redit$LLR(data = redit_data, groups = redit_groups)
                p <- get_numeric(llr, "p.value")
                if (is.na(p)) stop("REDIT-LLR did not return a valid p.value.")
                p
            })
            LLR_p <- llr_ans$p; LLR_error <- llr_ans$error

            reg_ans <- tryCatch({
                covariates <- data.frame(group = factor(redit_groups, levels = c(g1, g2)))
                reg <- redit$Regression(data = redit_data, covariates = covariates)
                effect_name <- paste0("group.", g2); p_name <- paste0(effect_name, ".p.value")
                p <- get_numeric(reg, p_name)
                if (is.na(p)) stop("REDIT-Regression did not return a valid group p.value.")
                list(p = p, effect = get_numeric(reg$parameter_estimates, effect_name),
                     intercept = get_numeric(reg$parameter_estimates, "intercept"),
                     sigma = get_numeric(reg$parameter_estimates, "sigma"),
                     converged = get_logical(reg$converged), error = NA_character_)
            }, error = function(e) list(p = NA_real_, effect = NA_real_, intercept = NA_real_,
                                        sigma = NA_real_, converged = NA, error = conditionMessage(e)))
            Regression_p <- reg_ans$p; Regression_effect <- reg_ans$effect
            Regression_intercept <- reg_ans$intercept; Regression_sigma <- reg_ans$sigma
            Regression_converged <- reg_ans$converged; Regression_error <- reg_ans$error
        }

        results[[i]] <- data.frame(
            ID = site, comparison = paste0(g1, "_vs_", g2), group1 = g1, group2 = g2,
            n_group1_valid = n1, n_group2_valid = n2,
            valid_samples_group1 = paste(valid_samples1, collapse = ";"),
            valid_samples_group2 = paste(valid_samples2, collapse = ";"),
            mean_group1 = mean1, mean_group2 = mean2, delta_mean_group2_minus_group1 = mean2 - mean1,
            pooled_group1 = pooled1, pooled_group2 = pooled2,
            delta_pooled_group2_minus_group1 = pooled2 - pooled1,
            ttest_p = ttest_p, fisher_p = fisher_p, LLR_p = LLR_p, Regression_p = Regression_p,
            Regression_effect = Regression_effect, Regression_intercept = Regression_intercept,
            Regression_sigma = Regression_sigma, Regression_converged = Regression_converged,
            ttest_status = test_status(ttest_p, ttest_error, ttest_eligible),
            fisher_status = test_status(fisher_p, fisher_error, fisher_eligible),
            LLR_status = test_status(LLR_p, LLR_error, redit_eligible),
            Regression_status = test_status(Regression_p, Regression_error, redit_eligible),
            ttest_error = ttest_error, fisher_error = fisher_error, LLR_error = LLR_error,
            Regression_error = Regression_error, stringsAsFactors = FALSE)
        if (i %% 1000 == 0 || i == n_sites) cat("Processed:", i, "/", n_sites, "\n")
    }

    # BH-FDR correction and significance calls
    out <- rbindlist(results, fill = TRUE)
    out$ttest_FDR <- bh_adjust(out$ttest_p); out$fisher_FDR <- bh_adjust(out$fisher_p)
    out$LLR_FDR <- bh_adjust(out$LLR_p); out$Regression_FDR <- bh_adjust(out$Regression_p)
    out$ttest_significant <- !is.na(out$ttest_FDR) & out$ttest_FDR < FDR_CUTOFF
    out$fisher_significant <- !is.na(out$fisher_FDR) & out$fisher_FDR < FDR_CUTOFF
    out$LLR_significant <- !is.na(out$LLR_FDR) & out$LLR_FDR < FDR_CUTOFF
    out$Regression_significant <- !is.na(out$Regression_FDR) & out$Regression_FDR < FDR_CUTOFF

    # Write result and summary tables
    comparison_name <- paste0(g1, "_vs_", g2); result_dir <- file.path(outdir, comparison_name)
    dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
    result_file <- file.path(result_dir, paste0(comparison_name, "_hypothesis_tests.tsv"))
    summary_file <- file.path(result_dir, paste0(comparison_name, "_summary.tsv"))
    fwrite(out, result_file, sep = "\t", na = "NA")
    summary <- data.frame(
        comparison = comparison_name, group1 = g1, group2 = g2,
        group1_samples = paste(samples1, collapse = ";"), group2_samples = paste(samples2, collapse = ";"),
        input_sites = nrow(out),
        sites_with_valid_samples_both_groups = sum(out$n_group1_valid > 0 & out$n_group2_valid > 0),
        ttest_tested = sum(out$ttest_status == "tested"), fisher_tested = sum(out$fisher_status == "tested"),
        LLR_tested = sum(out$LLR_status == "tested"), Regression_tested = sum(out$Regression_status == "tested"),
        ttest_significant_sites = sum(out$ttest_significant),
        fisher_significant_sites = sum(out$fisher_significant),
        LLR_significant_sites = sum(out$LLR_significant),
        Regression_significant_sites = sum(out$Regression_significant),
        min_coverage_inclusive = min_coverage, ttest_min_samples_per_group = TTTEST_MIN,
        fisher_min_samples_per_group = FISHER_MIN, redit_min_samples_per_group = REDIT_MIN,
        fdr_cutoff = FDR_CUTOFF)
    fwrite(summary, summary_file, sep = "\t", na = "NA")
    cat("\nComparison completed:", comparison_name, "\n"); print(summary)
    cat("\nResults:", result_file, "\nSummary:", summary_file, "\n")

    # Export both mean-level and pooled-level scatter PDFs
    plot_pair(out, g1, g2, result_dir, level_type = "mean")
    plot_pair(out, g1, g2, result_dir, level_type = "pooled")
    invisible(out)
}

# Main
main <- function() {
    args <- commandArgs(trailingOnly = TRUE); opt <- parse_args(args)
    FDR_CUTOFF <<- opt$fdr_cutoff
    dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(opt$outdir)) stop("Cannot create output directory: ", opt$outdir)
    opt$outdir <- normalizePath(opt$outdir, winslash = "/", mustWork = TRUE)
    cat("========================================\nRNA editing pairwise hypothesis testing\n========================================\n")
    cat("Input:", opt$input, "\nOutput:", opt$outdir,
        "\nCoverage threshold (inclusive):", opt$min_coverage,
        "\nBH-FDR cutoff:", FDR_CUTOFF, "\n")
    cat("\nDefined groups:\n")
    for (g in names(opt$groups)) cat(g, ": ", paste(opt$groups[[g]], collapse = ", "), "\n", sep = "")
    cat("\nRequested comparisons:\n")
    for (p in opt$comparisons) cat(p[1], "vs", p[2], "\n")
    counts <- read_counts(opt$input, opt$groups); redit <- load_redit()
    cat("\nStatistical settings:\nStudent's t-test minimum per group:", TTTEST_MIN,
        "\nFisher minimum per group:", FISHER_MIN, "\nREDIT minimum per group:", REDIT_MIN,
        "\nBH-FDR cutoff:", FDR_CUTOFF, "\nTop labels per panel:", TOP_LABELS, "\n")
    dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
    for (p in opt$comparisons) {
        run_pair(counts, p[1], p[2], opt$groups[[p[1]]], opt$groups[[p[2]]],
                 opt$min_coverage, redit, opt$outdir)
    }
    cat("\n========================================\nAll requested comparisons completed.\n========================================\n")
}

main()
