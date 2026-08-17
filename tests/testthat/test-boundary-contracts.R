context("AUDIT3 boundary contracts")

# Regression tests for the audit3 boundary conditions:
# A. q=0 entropy is independent of pseudocount (raw-support policy)
# B. non-natural log_base is rejected for multi-q spectra
# C. duplicated q within subject x condition is rejected deterministically
# D. TPM + effective-length double normalization is rejected

test_that("A: q=0 entropy uses raw support, independent of pseudocount", {
    x <- c(100, 0, 0)

    s0_no_pc <- .calculate_tsallis_entropy(x, q = 0, what = "S", pseudocount = 0)
    s0_pc <- .calculate_tsallis_entropy(x, q = 0, what = "S", pseudocount = 0.5)

    # Raw support is 1 observed isoform -> S_0 = 0, regardless of pseudocount
    expect_equal(unname(s0_pc), 0)
    expect_equal(unname(s0_pc), unname(s0_no_pc))

    # Hill D_0 = effective richness = raw support count = 1
    d0_pc <- .calculate_tsallis_entropy(x, q = 0, what = "D", pseudocount = 0.5)
    expect_equal(unname(d0_pc), 1)

    # .entropy_single follows the same policy
    e_no_pc <- .entropy_single(x, q = 0, norm = FALSE, pseudocount = 0)
    e_pc <- .entropy_single(x, q = 0, norm = FALSE, pseudocount = 0.5)
    expect_equal(e_pc, 0)
    expect_equal(e_pc, e_no_pc)
})

test_that("B: log_base != exp(1) is rejected for multi-q spectra", {
    x <- c(10, 5, 0)
    expect_error(
        .calculate_tsallis_entropy(x, q = c(1, 2), log_base = 2),
        "log_base"
    )
    expect_error(
        .calculate_tsallis_entropy(x, q = c(0.5, 1.5), log_base = 10),
        "log_base"
    )
    # single-q with non-natural base remains allowed (q=1 Shannon convention)
    expect_true(is.finite(.calculate_tsallis_entropy(x, q = 1, log_base = 2, norm = FALSE)))
})

test_that("C: duplicated q within subject x condition is rejected deterministically", {
    df <- data.frame(
        entropy = rnorm(8),
        q = c(0, 1, 2, 3, 0, 1, 2, 3),
        subject = factor(rep(1:2, each = 4)),
        condition = factor(rep(c("A", "B"), each = 4))
    )
    # duplicate q=3 within subject 1 x condition A block
    df_dup <- rbind(
        df,
        data.frame(
            entropy = 0,
            q = 3,
            subject = factor(1, levels = levels(df$subject)),
            condition = factor("A", levels = levels(df$condition))
        )
    )

    expect_error(
        .build_ar1_cor(df_dup, grid_col = "time_idx"),
        "Duplicated q values"
    )

    # clean data builds fine
    expect_true(is.list(.build_ar1_cor(df, grid_col = "time_idx")))
})

test_that("D: tpm = TRUE with effective_length is rejected (no double normalization)", {
    skip_if_not_installed("SummarizedExperiment")

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            counts = matrix(seq_len(12), nrow = 4),
            tpm = matrix(stats::runif(12), nrow = 4)
        ),
        colData = S4Vectors::DataFrame(sample = paste0("S", 1:3),
            condition = c("A", "A", "B"))
    )
    S4Vectors::metadata(se)$effective_length <- c(100, 200, 300, 400)
    rownames(se) <- paste0("TX", 1:4)
    genes <- c("G1", "G1", "G2", "G2")

    expect_error(
        .calculate_diversity(se, genes = genes, q = 1, tpm = TRUE, verbose = FALSE),
        "effective_length"
    )
})
