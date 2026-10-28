list("value$symmetric$internal" = list(covariate_names = c("grna_n_nonzero", 
"grna_n_umis", "response_n_nonzero", "response_n_umis", "response_p_mito"
), low_moi = TRUE, tables = list(grna_target_data_frame = list(
    class = "data.frame", names = c("grna_id", "grna_target"), 
    row_names = 1:2, columns = list(grna_id = c("g1", "g2"), 
        grna_target = c("non-targeting", "target"))))), "value$near_symmetric$internal" = list(
    covariate_names = c("grna_n_nonzero", "grna_n_umis", "response_n_nonzero", 
    "response_n_umis", "response_p_mito"), low_moi = TRUE, tables = list(
        grna_target_data_frame = list(class = "data.frame", names = c("grna_id", 
        "grna_target"), row_names = 1:2, columns = list(grna_id = c("g1", 
        "g2"), grna_target = c("non-targeting", "target"))))))
