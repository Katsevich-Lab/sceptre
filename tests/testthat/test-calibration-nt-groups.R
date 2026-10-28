test_that("calibration checks the group size against NT gRNAs remaining after QC", {
    inputs <- refactor_make_inputs()
    case <- refactor_case_definitions()$low_complement_perm
    object <- suppressMessages(refactor_prepare(inputs, case))
    groups <- object@grna_assignments_raw$indiv_nt_grna_idxs
    expect_true(length(groups) >= 2L)
    removed_group <- names(groups)[[1L]]
    # Remove every cell of one NT gRNA through QC, keeping its metadata row
    object <- suppressMessages(sceptre::run_qc(
        object,
        response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1),
        n_nonzero_trt_thresh = 0L,
        n_nonzero_cntrl_thresh = 0L,
        additional_cells_to_remove = unique(c(
            inputs$remove_cells,
            groups[[1L]]
        ))
    ))
    expect_false(
        removed_group %in% names(object@grna_assignments$indiv_nt_grna_idxs)
    )
    expect_true(removed_group %in% object@grna_target_data_frame$grna_id)
    active <- length(object@grna_assignments$indiv_nt_grna_idxs)
    expect_error(
        suppressMessages(sceptre::run_calibration_check(
            object,
            calibration_group_size = active + 1L,
            n_calibration_pairs = 1L
        )),
        "non-targeting gRNAs with cells remaining after QC"
    )
    result <- suppressMessages(sceptre::run_calibration_check(
        object,
        calibration_group_size = active,
        n_calibration_pairs = 1L,
        print_progress = FALSE
    ))
    expect_equal(nrow(result@calibration_result), 1L)

    one_active <- suppressMessages(sceptre::run_qc(
        object,
        response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1),
        n_nonzero_trt_thresh = 0L,
        n_nonzero_cntrl_thresh = 0L,
        additional_cells_to_remove = unique(c(
            inputs$remove_cells,
            unlist(groups[-length(groups)], use.names = FALSE)
        ))
    ))
    expect_length(one_active@grna_assignments$indiv_nt_grna_idxs, 1L)
    expect_gt(
        compute_calibration_group_size(one_active@grna_target_data_frame),
        1L
    )
    # The default group size comes from the metadata; check it against the
    # NT gRNAs that remain rather than silently reducing it
    expect_error(
        suppressMessages(sceptre::run_calibration_check(
            one_active,
            n_calibration_pairs = 1L
        )),
        "non-targeting gRNAs with cells remaining after QC"
    )
})

test_that("NT-cell calibration needs two NT gRNAs remaining after QC", {
    inputs <- refactor_make_inputs()
    case <- refactor_case_definitions()$low_nt_perm
    object <- suppressMessages(refactor_prepare(inputs, case))
    groups <- object@grna_assignments_raw$indiv_nt_grna_idxs
    one_active <- suppressMessages(sceptre::run_qc(
        object,
        response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1),
        n_nonzero_trt_thresh = 0L,
        n_nonzero_cntrl_thresh = 0L,
        additional_cells_to_remove = unique(c(
            inputs$remove_cells,
            unlist(groups[-length(groups)], use.names = FALSE)
        ))
    ))
    expect_length(one_active@grna_assignments$indiv_nt_grna_idxs, 1L)
    expect_error(
        suppressMessages(sceptre::run_calibration_check(
            one_active,
            calibration_group_size = 1L,
            n_calibration_pairs = 1L
        )),
        "Two or more non-targeting gRNAs with cells remaining after QC"
    )
})
