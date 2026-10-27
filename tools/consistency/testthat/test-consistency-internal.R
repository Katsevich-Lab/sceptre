# Internal state: sceptre_object slots (cached fits, budgets, QC'd index
# lists, the requested pairs with their pairwise QC). A restructuring that
# moves slots may waive these checks with `consistency.R run --waive-internal`,
# and must then say so in its pull request (README.md).
for (name in names(consistency_cases())) {
    local({
        case_name <- name
        test_that(paste("internal state:", case_name), {
            skip_if(
                consistency_internal_waived(),
                "internal-slot comparison waived (--waive-internal)"
            )
            consistency_expect_same(
                consistency_internal(consistency_actual(case_name)),
                consistency_internal(consistency_expected(case_name)),
                case_name,
                "Internal"
            )
        })
    })
}
