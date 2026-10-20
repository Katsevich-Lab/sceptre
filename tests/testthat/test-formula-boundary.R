test_that("generated formulas treat unusual covariate names as data labels", {
    names <- c("batch name", "if", "x+y", "tick`name", "+", "~", ".", "...", "..1", "offset", "offset dose")
    for (name in names) {
        covariates <- data.frame(value = factor(rep(c("a", "b"), 5L)),
            dropped = seq_len(10L), check.names = FALSE)
        names(covariates)[1L] <- name
        formula <- auto_construct_formula_object(covariates, FALSE)
        # dropped has only ten distinct values, so it remains a separate term;
        # a literal '.' column must not turn into an expansion of both terms.
        matrix <- stats::model.matrix(formula, covariates)
        expect_identical(dim(matrix), c(10L, 3L), info = name)
        expect_identical(unname(matrix[, 2L]), rep(c(0, 1), 5L), info = name)
        expect_identical(unname(matrix[, 3L]), as.numeric(seq_len(10L)), info = name)
        expect_identical(attr(matrix, "contrasts")[[1L]], "contr.treatment", info = name)
        expect_identical(environment(formula), baseenv())
    }
})

test_that("count transformations preserve log branches for literal unusual names", {
    covariates <- data.frame(a = c(0, 2, 3, 1), b = c(1, 4, 2, 5), check.names = FALSE)
    names(covariates) <- c("n_umis + special`name", "if n_nonzero")
    formula <- auto_construct_formula_object(covariates, TRUE)
    matrix <- stats::model.matrix(formula, covariates)
    expect_identical(unname(matrix[, 2L]), log(covariates[[1L]] + 1))
    expect_identical(unname(matrix[, 3L]), log(covariates[[2L]]))
    expect_identical(environment(formula), baseenv())
})

test_that("ordinary generated terms, design contrasts and RNG remain exact", {
    withr::local_seed(20261019)
    covariates <- data.frame(aaa_n_umis = rep(0:5, 4L),
        n_nonzero_aaa = seq_len(24L), batch = factor(rep(c("a", "b", "c"), 8L)),
        many = seq_len(24L), response_p_mito = rep(0.1, 24L),
        grna_n_umis = rep(c(4, 7, 9), 8L), grna_n_nonzero = rep(c(1, 2), 12L))
    before <- .Random.seed
    for (include in c(FALSE, TRUE)) {
        expected <- if (include) {
            ~ log(aaa_n_umis + 1) + log(n_nonzero_aaa) + batch + log(grna_n_umis) + log(grna_n_nonzero)
        } else ~ log(aaa_n_umis + 1) + log(n_nonzero_aaa) + batch
        actual <- auto_construct_formula_object(covariates, include)
        expect_identical(as.character(actual), as.character(expected))
        expect_identical(stats::model.matrix(actual, covariates), stats::model.matrix(expected, covariates))
        expect_identical(environment(actual), baseenv())
    }
    expect_identical(.Random.seed, before)
    no_terms <- data.frame(many = seq_len(20L))
    expect_identical(stats::model.matrix(auto_construct_formula_object(no_terms, FALSE), no_terms),
        stats::model.matrix(~ 1, no_terms))
})

test_that("user formula environments retain local transformations", {
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    case <- readRDS(test_path("fixtures", "refactor", "cases.rds"))$low_complement_perm
    object <- refactor_import(inputs, case)
    formula <- local({
        shifted <- function(x) x + 17
        ~ shifted(z) + batch
    })
    updated <- set_analysis_parameters(object, formula_object = formula,
        discovery_pairs = inputs$discovery_pairs)
    expect_identical(environment(updated@spec@formula_object), environment(formula))
    expect_identical(updated@preparation@design_matrix,
        stats::model.matrix(formula, object@data@covariate_data_frame))
})

test_that("offset labels remain usable while actual offset calls are refused", {
    expect_false(formula_contains_offset_call(~ offset + `offset dose`))
    expect_false(formula_contains_offset_call(~ log(`offset n_umis`)))
    expect_true(formula_contains_offset_call(~ offset(x)))
    expect_true(formula_contains_offset_call(~ stats::offset(x)))
    expect_true(formula_contains_offset_call(~ log(offset(x))))
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    case <- readRDS(test_path("fixtures", "refactor", "cases.rds"))$low_complement_perm
    object <- refactor_import(inputs, case)
    object@data@covariate_data_frame$offset <- sin(seq_len(nrow(object@data@covariate_data_frame)))
    expect_no_error(set_analysis_parameters(object, formula_object = ~ offset + batch))
    expect_error(set_analysis_parameters(object, formula_object = ~ offset(z)), "Offsets are not currently supported")
})


test_that("valid unusual imported names work through public default setup", {
    inputs <- readRDS(test_path("fixtures", "refactor", "inputs.rds"))
    for (name in c("batch name", "if", "tick`name", "+", "offset dose", ".", "...", "..1")) {
        covariates <- inputs$extra_covariates
        names(covariates)[names(covariates) == "batch"] <- name
        object <- import_data(response_matrix = inputs$response_matrix,
            grna_matrix = inputs$low_grna_matrix,
            grna_target_data_frame = inputs$grna_target_data_frame,
            extra_covariates = covariates, moi = "low")
        updated <- set_analysis_parameters(object, discovery_pairs = inputs$discovery_pairs)
        expect_identical(environment(updated@spec@formula_object), baseenv(), info = name)
        expect_true(stage_complete(updated, "design"), info = name)
        expect_identical(nrow(updated@preparation@design_matrix), ncol(inputs$response_matrix))
    }
})
