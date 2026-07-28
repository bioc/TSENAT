# ============================================================================
# TRANSCRIPT-LEVEL PLOTTING, VIOLIN/DENSITY & DIVERGENCE
# Extracted from plots_helpers.R — July 2026 refactoring (I11)
# ============================================================================


# ============================================================================
# PLOT EXPRESSION HELPERS
# ============================================================================
#' Make plot for a single gene
#'
#' Creates a transcript-level expression plot for a single gene.
#'
#' @param gene_single Character: gene identifier
#' @param mapping Data frame with Transcript and Gen columns
#' @param counts Matrix of read counts
#' @param samples Character vector of sample assignments
#' @param top_n Number of top transcripts to show
#' @param agg_fun Aggregation function for summarization
#' @param pseudocount Pseudocount for log transformation
#' @param agg_label_unique Label for aggregation metric
#' @param fill_limits Optional numeric vector for fill scale limits
#' @param font_scale Font scaling factor
#'
#' @return ggplot2 object
#' @noRd
.make_single_gene_plot <- function(gene_single, mapping, counts, samples,
    top_n, agg_fun, pseudocount, agg_label_unique, fill_limits = NULL, font_scale = 1) {
    built <- .build_transcript_long(gene_single, mapping, counts, samples,
        NULL)
    df_summary <- .aggregate_transcript_data(built$df_long, agg_fun, pseudocount)
    .build_tile_plot_from_summary(df_summary, agg_label_unique, fill_limits,
        font_scale = font_scale)
}

#' Combine multiple gene plots
#'
#' Combines individual gene plots into a grid layout.
#'
#' @param plots List of ggplot2 objects (one per gene)
#' @param output_file Optional file path to save combined plot
#' @param agg_label_unique Label for aggregation metric
#'
#' @return Combined plot object or invisible NULL if output_file provided
#' @noRd

.combine_gene_plots <- function(plots, output_file = NULL, agg_label_unique = NULL) {
    # Allow callers to pass a single character second argument as the
    # `agg_label_unique` for convenience (legacy test call patterns).
    if (is.null(agg_label_unique) && !is.null(output_file) && is.character(output_file) &&
        length(output_file) == 1) {
        agg_label_unique <- output_file
        output_file <- NULL
    }
    if (requireNamespace("patchwork", quietly = TRUE)) {
        .combine_plots_patchwork(plots, agg_label_unique)
    } else if (requireNamespace("cowplot", quietly = TRUE)) {
        .combine_plots_cowplot(plots, output_file = output_file, agg_label_unique = agg_label_unique)
    } else {
        .combine_plots_grid(plots, output_file = output_file, agg_label_unique = agg_label_unique)
    }
}


# ============================================================================
# TRANSCRIPT PLOTTING HELPERS: Data Preparation
# ============================================================================

#' Prepare Inputs for Transcript-Level Plotting
#'
#' Normalizes and validates counts, samples, and tx2gene mapping.
#' Creates aggregation function based on chosen metric.
#'
#' @param counts Matrix or data.frame with transcripts as rows, samples as
#' columns.
#'   Can also be a \code{SummarizedExperiment}.
#' @param readcounts Character: name of assay in SE (if counts is SE).
#' Default: NULL.
#' @param samples Character vector: sample group assignments (optional).
#' @param coldata Character/data.frame: sample metadata (optional).
#' @param condition_col Character: column name for grouping
#'   in coldata (default: 'sample_type').
#' @param tx2gene data.frame/character: Transcript-to-gene mapping with columns
#'   'Transcript' and 'Gen'. Can be file path or data.frame.
#' @param res Optional data.frame with results (gene names and p-values).
#' @param top_n Integer: number of transcripts to select.
#' @param pseudocount Numeric: pseudocount for log transformation (default: 0).
#' @param output_file Character: file path for saving plot (optional).
#' @param metric Character: aggregation metric
#'   ('median' [default], 'mean', 'variance', 'iqr').
#'
#' @return List with elements:
#'   - counts: normalized count matrix
#'   - samples: sample group assignments
#'   - mapping: tx2gene data.frame
#'   - metric_choice: chosen metric
#'   - agg_fun: aggregation function
#'   - agg_label_unique: metric label for display
#'   - top_n: number of transcripts
#'   - pseudocount: pseudocount value
#'   - output_file: output file path (if provided)
#'

#' @noRd

.prepare_transcript_inputs <- function(counts, readcounts = NULL, samples = NULL,
    coldata = NULL, condition_col = "condition", tx2gene = NULL, res = NULL, top_n = NULL,
    pseudocount = 0, output_file = NULL, metric = c("median", "mean", "variance",
        "iqr")) {


    # Handle SummarizedExperiment input
    if (inherits(counts, "SummarizedExperiment")) {
        se <- counts
        counts_mat <- .get_readcounts_from_se(se, readcounts)
        counts <- as.matrix(counts_mat)
        samples <- .infer_samples_from_se(se, samples, condition_col = condition_col)

        if (is.null(tx2gene)) {
            txres <- .get_tx2gene_from_se(se, counts)
            if (!is.null(txres) && !is.null(txres$mapping)) {
                mapping <- data.frame(Transcript = rownames(counts), Gen = as.character(txres$mapping),
                  stringsAsFactors = FALSE)
                tx2gene <- mapping
            }
        }
    }

    # Validate counts
    if (!is.matrix(counts) && !is.data.frame(counts)) {
        stop("`counts` must be a matrix, data.frame, or SummarizedExperiment", call. = FALSE)
    }

    counts <- as.matrix(counts)
    if (is.null(rownames(counts))) {
        stop("`counts` must have rownames (transcript identifiers)", call. = FALSE)
    }

    # Infer samples from coldata if needed
    if (is.null(samples)) {
        if (!is.null(coldata)) {
            samples <- .infer_samples_from_coldata(coldata, counts, condition_col)
        } else {
            stop("Either 'samples' or 'coldata' must be provided", call. = FALSE)
        }
    }

    # Validate and normalize tx2gene
    if (is.null(tx2gene)) {
        stop("`tx2gene` must be provided", call. = FALSE)
    }
    mapping <- .read_tx2gene(tx2gene)

    if (length(samples) != ncol(counts)) {
        stop("Length of `samples` must equal columns in `counts`", call. = FALSE)
    }

    # Create aggregation function
    metric_choice <- match.arg(metric)
    agg_fun <- switch(metric_choice, median = function(x) stats::median(x, na.rm = TRUE),
        mean = function(x) base::mean(x, na.rm = TRUE), variance = function(x) stats::var(x,
            na.rm = TRUE), iqr = function(x) stats::IQR(x, na.rm = TRUE))

    agg_label_metric <- if (metric_choice == "iqr")
        "IQR" else metric_choice
    agg_label <- sprintf("Transcript-level expression with metric %s", agg_label_metric)
    agg_label_unique <- agg_label

    list(counts = counts, samples = samples, mapping = mapping, metric_choice = metric_choice,
        agg_fun = agg_fun, agg_label_unique = agg_label_unique, top_n = top_n, pseudocount = pseudocount,
        output_file = output_file)
}

#' Read and Validate tx2gene Mapping
#'
#' Reads transcript-to-gene mapping from file or data.frame.
#' Validates required columns: 'Transcript' and 'Gen'.
#'
#' @param tx2gene Character (file path) or data.frame mapping.
#'
#' @return data.frame with columns 'Transcript' and 'Gen'.
#'

#' @noRd

.read_tx2gene <- function(tx2gene) {
    if (is.null(tx2gene)) {
        stop("`tx2gene` must be provided as file path or data.frame", call. = FALSE)
    }

    if (is.character(tx2gene) && length(tx2gene) == 1) {
        if (!file.exists(tx2gene)) {
            stop("tx2gene file not found: ", tx2gene, call. = FALSE)
        }
        mapping <- utils::read.delim(tx2gene, stringsAsFactors = FALSE, header = TRUE)
    } else if (is.data.frame(tx2gene)) {
        mapping <- tx2gene
    } else {
        stop("`tx2gene` must be file path or data.frame", call. = FALSE)
    }

    if (!all(c("Transcript", "Gen") %in% colnames(mapping))) {
        stop("tx2gene must have columns 'Transcript' and 'Gen'", call. = FALSE)
    }

    return(mapping)
}

#' Infer Samples from Column Metadata
#'
#' Extracts sample group assignments from coldata.
#' Aligns sample IDs from coldata to counts columns.
#'
#' @param coldata Character (file path) or data.frame with sample metadata.
#' @param counts Count matrix (for column name alignment).
#' @param condition_col Character: column name for grouping variable.
#'
#' @return Character vector of sample group assignments.
#'

#' @noRd

.infer_samples_from_coldata <- function(coldata, counts, condition_col) {
    if (is.character(coldata) && length(coldata) == 1) {
        if (!file.exists(coldata)) {
            stop("coldata file not found: ", coldata, call. = FALSE)
        }
        cdf <- utils::read.delim(coldata, header = TRUE, stringsAsFactors = FALSE)
    } else if (is.data.frame(coldata)) {
        cdf <- coldata
    } else {
        stop("`coldata` must be file path or data.frame", call. = FALSE)
    }

    # Try row-indexed matching first
    if (!is.null(rownames(cdf)) && all(colnames(counts) %in% rownames(cdf))) {
        return(as.character(cdf[colnames(counts), condition_col]))
    }

    # Try sample ID column matching
    sample_id_cols <- c("sample", "Sample", "sample_id", "id")
    sid <- intersect(sample_id_cols, colnames(cdf))

    if (length(sid) > 0) {
        sid <- sid[1]
        if (!all(colnames(counts) %in% as.character(cdf[[sid]]))) {
            stop("coldata sample ID column doesn't match counts columns", call. = FALSE)
        }
        row_ix <- match(colnames(counts), as.character(cdf[[sid]]))
        return(as.character(cdf[[condition_col]][row_ix]))
    }

    stop("Could not match coldata to counts. Provide rownames or sample ID column.",
        call. = FALSE)
}

#' Create Aggregation Function
#'
#' Builds an aggregation function based on chosen metric.
#'
#' @param metric Character: 'median' (default), 'mean', 'variance', or 'iqr'.
#'
#' @return List with:
#'   - metric_choice: the selected metric
#'   - agg_fun: function that computes the metric
#'   - agg_label_unique: display label
#'

#' @noRd

.create_aggregation_function <- function(metric = c("median", "mean", "variance",
    "iqr")) {
    metric_choice <- match.arg(metric)

    agg_fun <- switch(metric_choice, median = function(x) stats::median(x, na.rm = TRUE),
        mean = function(x) base::mean(x, na.rm = TRUE), variance = function(x) stats::var(x,
            na.rm = TRUE), iqr = function(x) stats::IQR(x, na.rm = TRUE))

    agg_label_metric <- if (metric_choice == "iqr")
        "IQR" else metric_choice
    agg_label <- sprintf("Transcript-level expression with metric %s", agg_label_metric)
    agg_label_unique <- agg_label

    list(metric_choice = metric_choice, agg_fun = agg_fun, agg_label_unique = agg_label_unique)
}

#' Build Long-Format Transcript Data
#'
#' Transforms wide count matrix to long-format data.frame
#' for ggplot visualization.
#'
#' @param gene_single Character: single gene identifier.
#' @param mapping data.frame: tx2gene mapping with Transcript and Gen columns.
#' @param counts Matrix: transcript count matrix.
#' @param samples Character vector: sample group assignments.
#' @param top_n Integer: limit to top N transcripts (optional).
#'
#' @return List with:
#'   - df_long: long-format data.frame (columns: tx, sample, expr, group)
#'   - txs: selected transcript identifiers
#'

#' @noRd

.build_transcript_long <- function(gene_single, mapping, counts, samples, top_n = NULL) {
    txs <- mapping$Transcript[mapping$Gen == gene_single]
    txs <- intersect(txs, rownames(counts))

    if (length(txs) == 0) {
        stop("No transcripts found for gene: ", gene_single, call. = FALSE)
    }

    if (!is.null(top_n)) {
        txs <- head(txs, top_n)
    }

    # Create long-format data
    mat <- counts[txs, , drop = FALSE]
    df_all <- as.data.frame(mat)
    df_all$tx <- rownames(mat)

    df_long <- tidyr::pivot_longer(df_all, -tx, names_to = "sample", values_to = "expr")
    df_long$group <- rep(samples, times = length(txs))

    list(df_long = df_long, txs = txs)
}

#' Aggregate Long-Format Transcript Data
#'
#' Summarizes expression by transcript and group.
#'
#' @param df_long Long-format data.frame from \code{build_transcript_long}.
#' @param agg_fun Function: aggregation function (e.g., median, mean).
#' @param pseudocount Numeric: pseudocount for log transformation.
#'
#' @return data.frame with columns: tx, group, expr, log2expr.
#'

#' @noRd

.aggregate_transcript_data <- function(df_long, agg_fun, pseudocount = 0) {
    df_summary <- stats::aggregate(expr ~ tx + group, data = df_long, FUN = agg_fun)
    df_summary$log2expr <- log2(df_summary$expr + pseudocount)
    df_summary$tx <- factor(df_summary$tx, levels = unique(df_summary$tx))
    return(df_summary)
}
# ============================================================================
# PLOT TSALLIS Q-CURVE HELPERS

# ============================================================================

#' Extract diversity objects from analysis results
#'
#' @param div_list List of diversity results (SummarizedExperiment or matrix)
#' @return List with: objects, q_names, first_se, bootstrap_ci_available

#' @noRd
.extract_diversity_objects <- function(div_list) {

    combined_assays_dict <- list()
    first_se <- NULL
    bootstrap_ci_available <- FALSE

    for (q_name in names(div_list)) {
        obj <- div_list[[q_name]]
        if (methods::is(obj, "SummarizedExperiment")) {
            if (is.null(first_se)) {
                first_se <- obj
            }
            # Check if bootstrap CIs are available
            si <- SummarizedExperiment::assayNames(obj)
            if ("ci_lower" %in% si && "ci_upper" %in% si) {
                bootstrap_ci_available <- TRUE
            }
        }
        mat <- if (methods::is(obj, "SummarizedExperiment")) {
            SummarizedExperiment::assay(obj, 1)
        } else {
            as.matrix(obj)
        }

        q_val <- as.numeric(sub("^q_", "", q_name))
        combined_assays_dict[[q_name]] <- list(matrix = mat, q_val = q_val, se_obj = obj)
    }

    if (is.null(first_se)) {
        stop("No valid SummarizedExperiment found in analysis@diversity_results")
    }

    list(objects = combined_assays_dict, q_names = names(combined_assays_dict), first_se = first_se,
        bootstrap_ci_available = bootstrap_ci_available)
}

#' Normalize matrix dimensions and row order
#'
#' @param matrix Matrix to normalize
#' @param target_genes Target gene order (character vector)
#' @param target_n_cols Target number of columns
#' @return Normalized matrix

#' @noRd
.normalize_matrix_to_target <- function(matrix, target_genes, target_n_cols) {
    # Adjust column count
    if (ncol(matrix) != target_n_cols) {
        if (ncol(matrix) > target_n_cols) {
            matrix <- matrix[, seq_len(target_n_cols), drop = FALSE]
        } else {
            pad_cols <- target_n_cols - ncol(matrix)
            matrix <- cbind(matrix, matrix(0, nrow = nrow(matrix), ncol = pad_cols))
        }
    }

    # Reorder rows to match target genes
    matrix[target_genes, , drop = FALSE]
}

#' Extract bootstrap CI matrices from SE object
#'
#' @param se_obj SummarizedExperiment or matrix object
#' @param target_genes Target gene order
#' @param target_n_cols Target number of columns
#' @param assay_names Assay names in SE
#' @return List(ci_lower, ci_upper) or NULL

#' @noRd
.extract_bootstrap_ci_matrices <- function(se_obj, target_genes, target_n_cols, assay_names) {

    if (!methods::is(se_obj, "SummarizedExperiment")) {
        return(NULL)
    }

    if (!("ci_lower" %in% assay_names && "ci_upper" %in% assay_names)) {
        return(NULL)
    }

    # Extract ci_lower
    ci_lower <- SummarizedExperiment::assay(se_obj, "ci_lower")
    ci_lower <- .normalize_matrix_to_target(ci_lower, target_genes, target_n_cols)

    # Extract ci_upper
    ci_upper <- SummarizedExperiment::assay(se_obj, "ci_upper")
    ci_upper <- .normalize_matrix_to_target(ci_upper, target_genes, target_n_cols)

    list(ci_lower = ci_lower, ci_upper = ci_upper)
}

#' Create Q-value suffixed column names
#'
#' @param colnames Column names (character vector or NULL)
#' @param q_val Q-value (numeric)
#' @param n_cols Number of column names needed
#' @return Character vector with _q=X.XXX suffix

#' @noRd
.create_q_suffixed_colnames <- function(colnames, q_val, n_cols) {
    if (is.null(colnames) || length(colnames) == 0) {
        colnames <- paste0("sample_", seq_len(n_cols))
    }

    clean_colnames <- sub("_q=.*$", "", colnames)
    paste0(clean_colnames, "_q=", formatC(q_val, format = "f", digits = 3))
}

#' Build combined colData across all q-values
#'
#' @param div_list Original diversity results list
#' @param q_names Q-value names (keys from div_list)
#' @param unique_colnames Final combined column names with q-suffix
#' @return Data frame with combined colData

#' @noRd
.build_combined_coldata <- function(div_list, q_names, unique_colnames_list) {

    combined_coldata_list <- list()

    for (q_name in q_names) {
        q_val <- as.numeric(sub("^q_", "", q_name))

        # Access unique_colnames by q-value name (stored as list keys in
        # .fill_combined_assays)
        unique_colnames <- unique_colnames_list[[q_name]]

        if (methods::is(div_list[[q_name]], "SummarizedExperiment")) {
            cd <- as.data.frame(SummarizedExperiment::colData(div_list[[q_name]]))
        } else {
            cd <- data.frame(row.names = unique_colnames)
        }

        cd$q <- q_val
        rownames(cd) <- unique_colnames
        combined_coldata_list[[q_name]] <- cd
    }

    do.call(rbind, combined_coldata_list)
}

#' Create combined SummarizedExperiment with assays and metadata
#'
#' @param combined_assay Main diversity assay matrix
#' @param combined_ci_lower CI lower matrix (optional)
#' @param combined_ci_upper CI upper matrix (optional)
#' @param combined_coldata ColData frame
#' @param first_se Template SE for rowData
#' @return SummarizedExperiment object

#' @noRd
.create_combined_se_object <- function(combined_assay, combined_ci_lower, combined_ci_upper,
    combined_coldata, first_se) {

    # Extract or create rowData, ensuring dimensions match combined_assay
    rd_combined <- tryCatch({
        rd_temp <- SummarizedExperiment::rowData(first_se)
        if (!is.null(rd_temp) && nrow(rd_temp) == nrow(combined_assay)) {
            # Ensure rownames match combined_assay
            rownames(rd_temp) <- rownames(combined_assay)
            rd_temp
        } else {
            NULL
        }
    }, error = function(e) NULL)

    if (is.null(rd_combined) || nrow(rd_combined) != nrow(combined_assay)) {
        rd_combined <- data.frame(gene_id = rownames(combined_assay), row.names = rownames(combined_assay),
            stringsAsFactors = FALSE)
    } else {
        # Ensure rownames match even if we're using extracted rowData
        rownames(rd_combined) <- rownames(combined_assay)
    }

    # Validate dimensions
    if (ncol(combined_assay) != nrow(combined_coldata)) {
        stop("Column mismatch: assay has ", ncol(combined_assay), " columns but colData has ",
            nrow(combined_coldata), " rows")
    }
    if (nrow(combined_assay) != nrow(rd_combined)) {
        stop("Row mismatch: assay has ", nrow(combined_assay), " rows but rowData has ",
            nrow(rd_combined), " rows")
    }

    # Validate names match
    if (!identical(colnames(combined_assay), rownames(combined_coldata))) {
        stop("Column name mismatch between assay and colData")
    }
    if (!identical(rownames(combined_assay), rownames(rd_combined))) {
        stop("Row name mismatch between assay and rowData")
    }

    # Build assays list
    assays_list <- list(diversity = combined_assay)

    if (!is.null(combined_ci_lower) && !is.null(combined_ci_upper)) {
        ci_lower_valid <- sum(!is.na(combined_ci_lower)) > 0
        ci_upper_valid <- sum(!is.na(combined_ci_upper)) > 0

        if (ci_lower_valid && ci_upper_valid) {
            assays_list$ci_lower <- combined_ci_lower
            assays_list$ci_upper <- combined_ci_upper
        }
    }

    # Create SE
    combined_se <- SummarizedExperiment::SummarizedExperiment(assays = assays_list,
        colData = combined_coldata, rowData = rd_combined)

    # Add metadata if CI available
    if (!is.null(combined_ci_lower)) {
        S4Vectors::metadata(combined_se)$bootstrap_ci_count <- sum(!is.na(combined_ci_lower))
        S4Vectors::metadata(combined_se)$has_bootstrap_ci <- (sum(!is.na(combined_ci_lower)) >
            0)
    }

    combined_se
}

#' Prepare single q-value data for combined assays
#'
#' @param q_name Q-value name (key from combined_assays_dict)
#' @param combined_assays_dict Dictionary of matrices and metadata
#' @param target_genes Target gene order
#' @param target_n_cols Target columns per q-value
#' @param bootstrap_ci_available Boolean: CIs available
#' @return List with: unique_colnames, ncol, ci_lower, ci_upper

#' @noRd
.prepare_q_value_for_combining <- function(q_name, combined_assays_dict, target_genes,
    target_n_cols, bootstrap_ci_available) {

    mat <- combined_assays_dict[[q_name]]$matrix
    q_val <- combined_assays_dict[[q_name]]$q_val
    se_obj <- combined_assays_dict[[q_name]]$se_obj

    # Normalize matrix dimensions
    mat <- .normalize_matrix_to_target(mat, target_genes, target_n_cols)

    # Create q-suffixed column names
    unique_colnames <- .create_q_suffixed_colnames(colnames(mat), q_val, ncol(mat))

    # Extract CI matrices if available
    ci_lower <- ci_upper <- NULL
    if (bootstrap_ci_available) {
        sim_names <- SummarizedExperiment::assayNames(se_obj)
        ci_matrices <- .extract_bootstrap_ci_matrices(se_obj, target_genes, target_n_cols,
            sim_names)
        if (!is.null(ci_matrices)) {
            ci_lower <- ci_matrices$ci_lower
            ci_upper <- ci_matrices$ci_upper
        } else {
            ci_lower <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
            ci_upper <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
        }
    }

    list(matrix = mat, unique_colnames = unique_colnames, ncol_val = ncol(mat), ci_lower = ci_lower,
        ci_upper = ci_upper)
}

#' Fill combined assay matrices with data from all q-values
#'
#' @param combined_assays_dict Dictionary of matrices and metadata per q-value
#' @param q_names Q-value names (keys)
#' @param target_genes Target gene order
#' @param target_n_cols Target number of columns
#' @param bootstrap_ci_available Boolean: CIs available
#' @return List with: combined_assay, combined_ci_lower, combined_ci_upper, unique_colnames_list

#' @noRd
.fill_combined_assays <- function(combined_assays_dict, q_names, target_genes, target_n_cols,
    bootstrap_ci_available) {
    total_cols <- target_n_cols * length(q_names)

    # Initialize matrices
    combined_assay <- matrix(0, nrow = length(target_genes), ncol = total_cols)
    rownames(combined_assay) <- target_genes

    combined_ci_lower <- if (bootstrap_ci_available) {
        matrix(NA, nrow = length(target_genes), ncol = total_cols)
    } else NULL
    combined_ci_upper <- if (bootstrap_ci_available) {
        matrix(NA, nrow = length(target_genes), ncol = total_cols)
    } else NULL

    unique_colnames_list <- list()
    col_idx <- 1

    for (q_name in q_names) {
        result <- .prepare_q_value_for_combining(q_name, combined_assays_dict, target_genes,
            target_n_cols, bootstrap_ci_available)

        ncol_q <- result$ncol_val
        if (col_idx + ncol_q - 1 > total_cols) {
            stop("Dimension mismatch: ", col_idx, " to ", col_idx + ncol_q - 1, " exceeds total_cols=",
                total_cols)
        }

        # Fill main assay
        combined_assay[, col_idx:(col_idx + ncol_q - 1)] <- result$matrix
        unique_colnames_list[[q_name]] <- result$unique_colnames

        # Fill CI matrices if available
        if (bootstrap_ci_available && !is.null(result$ci_lower)) {
            combined_ci_lower[, col_idx:(col_idx + ncol_q - 1)] <- result$ci_lower
            combined_ci_upper[, col_idx:(col_idx + ncol_q - 1)] <- result$ci_upper
        }

        col_idx <- col_idx + ncol_q
    }

    list(combined_assay = combined_assay, combined_ci_lower = combined_ci_lower,
        combined_ci_upper = combined_ci_upper, unique_colnames_list = unique_colnames_list)
}

#' Convert TSENATAnalysis to combined SummarizedExperiment
#'
#' @param analysis TSENATAnalysis object with diversity_results
#' @return SummarizedExperiment with combined assay across all q-values

#' @noRd
.prepare_combined_se <- function(analysis) {

    div_list <- analysis@diversity_results

    # Step 1: Extract diversity objects and metadata
    extracted <- .extract_diversity_objects(div_list)

    # Step 2: Get target dimensions
    target_genes <- rownames(extracted$first_se)
    target_n_cols <- ncol(extracted$first_se)

    # Step 3: Fill combined assays
    filled <- .fill_combined_assays(extracted$objects, extracted$q_names, target_genes,
        target_n_cols, extracted$bootstrap_ci_available)

    # Step 4: Build combined colData (which defines the sample names via
    # rownames)
    combined_coldata_df <- .build_combined_coldata(div_list, extracted$q_names, filled$unique_colnames_list)

    # Step 5: Set column names on all assays to match colData rownames
    combined_colnames <- rownames(combined_coldata_df)
    colnames(filled$combined_assay) <- combined_colnames
    if (!is.null(filled$combined_ci_lower) && !is.null(filled$combined_ci_upper)) {
        colnames(filled$combined_ci_lower) <- combined_colnames
        colnames(filled$combined_ci_upper) <- combined_colnames
    }

    # Step 6: Create and return combined SE
    .create_combined_se_object(filled$combined_assay, filled$combined_ci_lower, filled$combined_ci_upper,
        combined_coldata_df, extracted$first_se)
}

#' Compute gene-level statistics (median +/- SD) by group and q-value
#'
#' @param long_data Long-format data frame with Gene, q, group, tsallis columns
#' @return Data frame with central tendency and spread by gene, group, q

#' @noRd
.compute_gene_group_stats <- function(long_data, metric = "iqr") {

    metric <- match.arg(tolower(metric), c("iqr", "sd"))
    long_data$qnum <- as.numeric(as.character(long_data$q))

    # Calculate spread based on metric choice
    if (metric == "iqr") {
        spread_calc <- quote(stats::IQR(tsallis, na.rm = TRUE)/2)
    } else {
        spread_calc <- quote(sqrt(stats::var(tsallis, na.rm = TRUE)))
    }

    dplyr::summarise(dplyr::group_by(long_data, group, qnum), central = median(tsallis,
        na.rm = TRUE), spread = !!spread_calc, .groups = "drop")
}

#' Volcano plot for differential results
#'
#' Create a volcano plot showing fold-change (x-axis) versus adjusted
#' p-value significance (y-axis). The function auto-detects a suitable x-axis
#' column if one is not provided and expects an adjusted p-value column for
#' significance coloring.
#'
#' Combine Volcano and MA-Tsallis Plots in a Grid Layout
#'
#' Creates a side-by-side grid layout with a volcano plot on the left and an
#' MA-Tsallis plot on the right.
#' Both plots are generated from differential analysis results data.
#'
#' @param diff_df Data.frame from differential analysis containing required
#' columns for both volcano and MA plots.
#' @param x_col Column name for x-axis in volcano plot (e.g.,
#' 'mean_difference'). Auto-detected if NULL.
#' @param padj_col Column name for adjusted p-values (default: 'padj').
#' @param label_thresh Threshold for volcano plot labels (default: 0.1).
#' @param sig_alpha Numeric significance threshold for adjusted p-values
#' (default: 0.05).
#' @param top_n Number of top genes to annotate in volcano plot (default: 5).
#' @param title_volcano Title for volcano plot. If NULL, auto-generated.
#' @param title_ma Title for MA plot (default: 'Tsallis-based MA plot').
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return A `ggplot2` object showing a 1x2 grid with volcano plot on the
#' left and MA plot on the right.
#'
#' @examples
#' # Simulate differential analysis results
#' x <- data.frame(
#'   genes = paste0('g', seq_len(20)),
#'   mean_difference = rnorm(20, sd = 1),
#'   padj = runif(20, 1e-5, 0.1),
#'   log2_fold_change = rnorm(20, sd = 0.8)
#' )



#' Internal helper to compute fill limits across multiple genes (not exported)
#' @noRd
.plot_transcript_fill_limits <- function(genes, mapping, counts, samples, top_n,
    agg_fun, pseudocount) {
    mins <- maxs <- c()
    for (g in genes) {
        txs <- mapping$Transcript[mapping$Gen == g]
        txs <- intersect(txs, rownames(counts))
        if (length(txs) == 0)
            next
        if (!is.null(top_n))
            txs <- head(txs, top_n)
        mat <- counts[txs, , drop = FALSE]
        df_all <- as.data.frame(mat)
        df_all$tx <- rownames(mat)
        df_long <- tidyr::pivot_longer(df_all, -tx, names_to = "sample", values_to = "expr")
        df_long$group <- rep(samples, times = length(txs))
        df_summary <- stats::aggregate(expr ~ tx + group, data = df_long, FUN = agg_fun)
        df_summary$log2expr <- log2(df_summary$expr + pseudocount)
        mins <- c(mins, min(df_summary$log2expr, na.rm = TRUE))
        maxs <- c(maxs, max(df_summary$log2expr, na.rm = TRUE))
    }
    if (length(mins) == 0)
        stop("No transcripts found for provided genes")
    c(min(mins, na.rm = TRUE), max(maxs, na.rm = TRUE))
}

#' Internal helper to draw grid layout with title, plots, and legend using
#' base grid
#' @noRd
.plot_transcript_grid_draw <- function(grobs, title, legend_grob, ncol, heights,
    to_file = NULL) {
    # If no output file is provided and no graphics device is open, render to a
    # temporary pdf device so that plotting in non-interactive sessions does
    # not create `Rplots.pdf` in the working directory.
    temp_dev <- FALSE
    # Only open a temporary PDF device when: - caller did not supply an output
    # file (`to_file` is NULL), - the session is non-interactive, and - no
    # graphics device is currently open (dev.cur() == 1)
    if (is.null(to_file) && !interactive() && grDevices::dev.cur() == 1L) {
        tmp <- tempfile("TSENAT_plot_", fileext = ".pdf")
        grDevices::pdf(tmp)
        temp_dev <- TRUE
        # Ensure device is closed and temporary file removed on exit
        on.exit({
            try(grDevices::dev.off(), silent = TRUE)
            if (file.exists(tmp)) unlink(tmp)
        }, add = TRUE)
    }

    # Calculate number of rows needed for plots
    nrow_plots <- ceiling(length(grobs)/ncol)
    nrow_total <- 2 + nrow_plots  # title + plot rows + legend

    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow_total, ncol,
        heights = heights)))
    # Title row
    vp_title <- grid::viewport(layout.pos.row = 1, layout.pos.col = seq_len(ncol))
    grid::pushViewport(vp_title)
    grid::grid.text("Transcript level expression", x = 0.5, y = 0.6, gp = grid::gpar(fontsize = 14,
        fontface = "bold"))
    grid::grid.text(paste0("Top genes with metric ", title), x = 0.5, y = 0.2, gp = grid::gpar(fontsize = 11,
        fontface = "italic", col = "gray40"))
    grid::upViewport()
    # Plot rows
    for (i in seq_along(grobs)) {
        plot_row_idx <- ((i - 1)%/%ncol) + 2
        plot_col_idx <- ((i - 1)%%ncol) + 1
        vp <- grid::viewport(layout.pos.row = plot_row_idx, layout.pos.col = plot_col_idx)
        grid::pushViewport(vp)
        grid::grid.draw(grobs[[i]])
        grid::upViewport()
    }
    # Legend row
    if (!is.null(legend_grob)) {
        vp_leg <- grid::viewport(layout.pos.row = nrow_total, layout.pos.col = seq_len(ncol))
        grid::pushViewport(vp_leg)
        grid::grid.draw(legend_grob)
        grid::upViewport()
    }
    grid::upViewport()
    # If caller supplied a file (caller likely opened a device), close it here.
    if (!is.null(to_file))
        grDevices::dev.off()
    invisible(NULL)
}

## Internal helpers for `plot_top_transcripts` refactor Create per-gene plot
## and combine multiple gene plots into final output
## NOTE: Primary definitions are at top of this file (lines ~146-166).

## Prepare and validate inputs for `plot_top_transcripts`

.build_tile_plot_from_summary <- function(df_summary, agg_label_unique,
    fill_limits = NULL, font_scale = 1) {
    # Calculate font sizes proportionally to output dimensions Reference: 12x8
    # inches (96 sq in) uses font_base=11 For other sizes, scale base font as:
    # base_font = 11 * sqrt(area/96) This ensures readability is maintained
    # across different output sizes

    base_font <- 11 * font_scale
    y_axis_font <- 12 * font_scale
    x_axis_font <- 14 * font_scale
    title_font <- 16 * font_scale
    legend_font <- 9 * font_scale

    p <- ggplot2::ggplot(df_summary, ggplot2::aes(x = group, y = tx, fill = log2expr)) +
        ggplot2::geom_tile(color = "black", linewidth = 0.3, width = 0.95, height = 0.92) +
        ggplot2::geom_vline(xintercept = 1.5, color = "white", linewidth = 1.5) +
        ggplot2::scale_x_discrete(expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = c(0,
        0)) + ggplot2::scale_fill_distiller(palette = "Blues", na.value = "lightgray",
        limits = fill_limits, name = "log2(expr)") + .theme_base(base_size = base_font) +
        ggplot2::labs(title = agg_label_unique, x = NULL, y = NULL, fill = "log2(expr)") +
        ggplot2::theme(plot.title = ggplot2::element_text(size = title_font, hjust = 0.5,
            face = "bold"), plot.margin = ggplot2::margin(4, 4, 4, 4))

    # Apply axis label formatting and legend configuration using Phase 5
    # helpers
    p <- .format_axis_labels(p, x_size = x_axis_font, y_size = y_axis_font, y_face = "plain",
        bold_title = FALSE)
    p <- .configure_legend(p, position = "bottom", width_cm = 2, text_size = legend_font)
    p <- p + ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top",
        barwidth = 10, barheight = 0.5, title.theme = ggplot2::element_text(size = title_font)))
    p
}
#' @importFrom ggplot2 ggplot aes geom_col geom_point geom_line scale_y_continuous labs theme_minimal theme element_text geom_hline geom_vline scale_color_manual scale_shape_manual geom_tile scale_fill_gradient2
NULL


# Internal plot helpers

.format_label <- function(lbl) {
    if (is.null(lbl)) {
        return(NULL)
    }
    s <- gsub("_", " ", lbl)
    s <- gsub("\\s+", " ", s)
    s <- trimws(s)
    s <- tolower(s)
    if (nchar(s) == 0) {
        return(s)
    }
    if (nchar(s) == 1) {
        return(toupper(s))
    }
    paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
}

# [REMOVED] .prepare_ma_plot_df() - MA plot functionality removed (April 2026)

#' Plot Tsallis Divergence Effect Size Distribution
#'
#' Generate a histogram visualization of Tsallis divergence effect sizes
#' across genes,
#' showing the distribution of information-theoretic measures of isoform
#' switching.
#'
#' @param interaction_results A data frame containing LMM results merged
#' with per-q divergence estimates.
#' Must contain columns matching the pattern `effect_size_D_q*` (e.g.,
#' `effect_size_D_q0_5`, `effect_size_D_q1_0`).
#'   Typically the result from [.calculate_effect_sizes()].
#'
#' @param threshold Numeric. Effect size threshold for visual marking.
#' Default is 0.1 (information-theoretic significance level).
#'
#' @return If ggplot2 is available and `interaction_results` contains valid
#' data, returns a ggplot object.
#'   Otherwise returns NULL invisibly and prints an informative message.
#'
#' @details
#' The function visualizes the distribution of effect sizes using the median
#' q-value's divergence
#' (typically around q=1.0, close to Shannon entropy). The red dashed line
#' marks the default
#' information-theoretic significance threshold of D=0.1.
#'
#' @references
#' - Chanda et al. (2020). Information Theory in Computational Biology.
#' *Entropy*, 22(6), 627.
#' - Tsallis, C. (1988). Possible Generalization of Boltzmann-Gibbs
#' Statistics. *Journal of Statistical Physics*, 52(1), 479-487.
#'
#' @examples
#' # Create example interaction results with divergence effect sizes
#' set.seed(123)
#' interaction_results <- data.frame(
#'   gene = paste0('gene_', 1:20),
#'   effect_size_D_q0_5 = runif(20, 0, 0.3),
#'   effect_size_D_q1_0 = runif(20, 0, 0.2),
#'   effect_size_D_q1_5 = runif(20, 0, 0.25)
#' )
#' 
#' # Plot divergence distribution
#' .plot_divergence_distribution(interaction_results, threshold = 0.1)
#'

#' @noRd

.plot_divergence_distribution <- function(interaction_results, threshold = 0.1) {

    # Check for ggplot2
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        message("ggplot2 package required for plotting. Please install: install.packages('ggplot2')")
        return(invisible(NULL))
    }

    # Validate input
    if (is.null(interaction_results) || nrow(interaction_results) == 0) {
        message("Plot not generated: interaction_results is empty or NULL.")
        return(invisible(NULL))
    }

    # Find per-q effect size columns
    effect_cols <- grep("^effect_size_D_q", colnames(interaction_results), value = TRUE)

    if (length(effect_cols) == 0) {
        message("Plot not generated: no per-q effect size columns found in interaction_results.")
        message("Expected columns like 'effect_size_D_q0_5', 'effect_size_D_q1_0', etc.")
        return(invisible(NULL))
    }

    # Use median q effect size for visualization
    median_idx <- ceiling(length(effect_cols)/2)
    median_col <- effect_cols[median_idx]

    # Filter out NA/NaN/Inf values for plotting
    plot_data <- interaction_results[is.finite(interaction_results[[median_col]]),
        , drop = FALSE]

    if (nrow(plot_data) == 0) {
        message("Plot not generated: no valid (finite) effect sizes to plot.")
        return(invisible(NULL))
    }

    # Create visualization of effect size distribution with publication theme
    p_effect <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[median_col]])) +
        ggplot2::geom_histogram(binwidth = 0.02, fill = .palette_blue_red()[1], alpha = 0.7,
            color = "black") + ggplot2::labs(title = "Distribution of Tsallis Divergence (D[q]) effect sizes across genes",
        subtitle = "Information-theoretic measure respecting Tsallis multi-q entropy properties",
        x = bquote("Effect size (Tsallis Divergence" ~ D[q] ~ "; D >" ~ .(threshold) ~
            "= meaningful information separation)"), y = "Number of genes", caption = paste("Red dashed line: D =",
            threshold, "filtering threshold (information-theoretic significance for q-dependent entropy)")) +
        .theme_base(base_size = 12) + ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "gray90"))

    # Add reference line using Phase 5 helper
    p_effect <- .add_reference_lines(p_effect, v_intercept = threshold, v_color = "red",
        v_size = 1)

    # Add threshold annotation
    p_effect <- p_effect + ggplot2::annotate("text", x = threshold, y = Inf, label = paste("Information\nthreshold\n(D=",
        threshold, ")", sep = ""), vjust = 1.5, hjust = 0, color = "red", size = 3.5)

    return(p_effect)
}




