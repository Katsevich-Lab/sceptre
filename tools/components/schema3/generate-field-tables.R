# Run from the package root after installing developer dependencies. Reflection
# is authoritative; missing or stale descriptions stop generation.
pkgload::load_all(".", quiet = TRUE, recompile = FALSE)
classes <- c("sceptre_object", "ScreenData", "AnalysisSpec", "PreparationRecord",
    "FitCache", "ResultRecord", "RunPlan", "CellView", "BackedCountStore")
metadata <- read.delim("tools/components/schema3/field-descriptions.tsv", check.names = FALSE)
output <- c("# Final schema 3 object fields", "",
    "Generated from `methods::getSlots()` by `generate-field-tables.R`. Declared",
    "types are exact. These tables describe the live schema, not legacy aliases.", "",
    "`CountStore` is the union of `dgRMatrix`, `BackedCountStore` and `NULL`.",
    "`CountHandle` is the union of `externalptr` and `NULL`. Neither union names",
    "an optional dependency class; loading the package does not load ondisc.", "")
actual <- list()
for (class in classes) {
    slots <- methods::getSlots(class)
    descriptions <- metadata[metadata$class == class, ]
    stopifnot(setequal(names(slots), descriptions$field), !anyDuplicated(descriptions$field))
    rows <- data.frame(class = class, field = names(slots), declared_type = unname(slots),
        description = descriptions$description[match(names(slots), descriptions$field)])
    actual[[class]] <- rows
    output <- c(output, paste0("## ", class), "",
        "| Field | Exact declared type | Meaning |", "|---|---|---|",
        paste0("| `", rows$field, "` | `", rows$declared_type, "` | ", rows$description, " |"), "")
}
writeLines(output, "tools/components/schema3/OBJECT_FIELDS.md")
write.table(do.call(rbind, actual), "tools/components/schema3/current-fields.tsv",
    sep = "\t", quote = FALSE, row.names = FALSE)
