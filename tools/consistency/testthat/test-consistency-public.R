# Public outputs: exported getters, get_result(), printed summaries, written
# outputs, pair constructors, errors, warnings and the RNG state. These must
# match for every change that is not meant to alter results.
for (name in names(consistency_cases())) {
    local({
        case_name <- name
        test_that(paste("public outputs:", case_name), {
            consistency_expect_same(
                consistency_public(consistency_actual(case_name)),
                consistency_public(consistency_expected(case_name)),
                case_name,
                "Public"
            )
        })
    })
}
