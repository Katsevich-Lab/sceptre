# A result boundary owns the entire exposed graph, including nested S4 slots
# and data.table attributes. Serialization preserves the exact ALTREP, compact
# row-name and attribute representations used by existing RunPlan certificates.
# Leave data.table's serialized null self-reference intact: its normal readRDS
# path allocates locally on the first := write, without a shallow-copy warning.
# Eager alloc.col repair would change serialized growable-vector flags on R4.4.
owned_result_value <- function(value) unserialize(serialize(value, NULL))
