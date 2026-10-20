from pathlib import Path
import csv,json
b=Path('/private/tmp/sceptre-submission-matched-check')
a=list(csv.DictReader(open(b/'submitted-functions.csv')));c=list(csv.DictReader(open(b/'qualified-functions.csv')))
a_idx={(x['file'],x['name']):x for x in a}
renamed={
('R/analysis.R','analyze_association_pairs'):('R/s4_analysis_functs_2.R','run_sceptre_analysis_high_level','0b6ae38 and f4aaab4: same public analysis dispatcher, permutation draw bank and eligible/failed-pair assembly'),
('R/association_execution.R','run_permutation_analysis'):('R/medium_level_functs_v2.R','run_perm_test_in_memory','0b6ae38: explicit old/new function replacement; context replaces argument expansion'),
('R/association_execution.R','analyze_response_block'):('R/medium_level_functs_v2.R','analyze_given_response_ids','0b6ae38: nested per-response block replacement; same response fitting and group testing loop'),
('R/association_execution.R','run_crt_analysis'):('R/medium_level_functs_v2.R','run_crt_in_memory_v2','0b6ae38 and ec1980c: explicit old/new CRT coordinator replacement; later compact perturbation fits added'),
('R/association_execution.R','analyze_guide_block'):('R/medium_level_functs_v2.R','perform_association_analysis','0b6ae38: explicit nested group-analysis replacement; guide groups retain request order'),
('R/association_workers.R','test_crt_shared_fit'):('R/association_testing_functions_v2.R','crt_glm_factored_out','0b6ae38: explicit replacement of shared response-fit CRT test; same per-guide sampling then per-response native tests'),
('R/association_workers.R','test_crt_group_fit'):('R/association_testing_functions_v2.R','discovery_ntcells_crt','0b6ae38: explicit replacement of treatment-plus-NT CRT test'),
('R/neg_control_functions.R','construct_negative_control_pairs'):('R/neg_control_functions.R','construct_negative_control_pairs_v2','Retained negative-control pair construction, suffix removal and storage/count adapters'),
('R/result_processing.R','pair_results_from_native'):('R/aux_functions.R','construct_data_frame_v2','f4aaab4: raw result field extraction moved verbatim then alias removed'),
('R/result_processing.R','integrate_guides'):('R/s4_analysis_functs_2.R','apply_grouping_to_result','f4aaab4: guide integration, metadata join and Bonferroni body extracted; adapter alias later removed'),
('R/plotting_functions.R','_anonymous_.157'):('R/plotting_functions.R','_anonymous_.154','Same anonymous lapply(grnas_to_plot, function(curr_grna_id)) graph-construction site'),
('R/s4_helpers.R','_anonymous_.87'):('R/s4_helpers.R','_anonymous_.83','Same setMethod(print, sceptre_object, function(x)); component slot and results access adapted'),
('R/storage.R','set_matrix_accessibility'):('R/preprocessing_functions.R','set_matrix_accessibility','0acfe2a: explicit same-name relocation into storage module'),
('R/s4_helpers.R','_anonymous_.2'):('R/s4_helpers.R','_anonymous_.2','Same setMethod(show, sceptre_object, function(object)); previously exactly 50 lines'),
}
rows=[]
for x in c:
 if int(x['length'])<=50:continue
 key=(x['file'],x['name']);r=renamed.get(key)
 if r:
  y=a_idx[r[:2]];basis=r[2]
 else:
  y=a_idx.get(key);basis='Same named function in same source file' if y else 'No one-to-one pre-existing long-function site; extracted or new contract/helper'
 if y and int(y['length'])>50:status='inherited_function_site'
 elif y:status='new_threshold_crossing'
 else:status='new_function_site'
 rows.append(dict(x,baseline_file=y['file'] if y else '',baseline_name=y['name'] if y else '',baseline_length=y['length'] if y else '',baseline_start_line=y['startLine'] if y else '',baseline_end_line=y['endLine'] if y else '',status=status,basis=basis))
with open(b/'qualified-long-function-lineage.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
from collections import Counter
print(Counter(x['status'] for x in rows))
for x in rows:
 if x['status']!='inherited_function_site':print(x['file']+':'+x['startLine'],x['name'],x['length'],x['status'],'baseline',x['baseline_length'])
