"""Conservative instance gate against the actual submitted package source.
Args: submitted extracted root, candidate extracted root, candidate-functions CSV,
      old submitted-functions CSV, output directory.
"""
import csv,json,sys,collections
from pathlib import Path
submitted,candidate,functions,baseline,out=map(Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=True)
read=lambda p:list(csv.DictReader(p.open()))
old=read(baseline);new=read(functions)
renames={
 ('R/analysis.R','analyze_association_pairs'):('R/s4_analysis_functs_2.R','run_sceptre_analysis_high_level'),
 ('R/association_execution.R','run_permutation_analysis'):('R/medium_level_functs_v2.R','run_perm_test_in_memory'),
 ('R/association_execution.R','analyze_response_block'):('R/medium_level_functs_v2.R','analyze_given_response_ids'),
 ('R/association_execution.R','run_crt_analysis'):('R/medium_level_functs_v2.R','run_crt_in_memory_v2'),
 ('R/association_execution.R','analyze_guide_block'):('R/medium_level_functs_v2.R','perform_association_analysis'),
 ('R/association_workers.R','test_crt_shared_fit'):('R/association_testing_functions_v2.R','crt_glm_factored_out'),
 ('R/association_workers.R','test_crt_group_fit'):('R/association_testing_functions_v2.R','discovery_ntcells_crt'),
 ('R/neg_control_functions.R','construct_negative_control_pairs'):('R/neg_control_functions.R','construct_negative_control_pairs_v2'),
 ('R/result_processing.R','pair_results_from_native'):('R/aux_functions.R','construct_data_frame_v2'),
 ('R/result_processing.R','integrate_guides'):('R/s4_analysis_functs_2.R','apply_grouping_to_result'),
 ('R/storage.R','set_matrix_accessibility'):('R/preprocessing_functions.R','set_matrix_accessibility'),
}
lineage=[]
for item in new:
 if int(item['length'])<=50:continue
 key=(item['file'],item['name']);target=renames.get(key,key);basis='same named site' if target==key else 'reviewed rename or relocation'
 if item['name'].startswith('_anonymous_'):
  first=(candidate/item['file']).read_text().splitlines()[int(item['startLine'])-1]
  if item['file']=='R/s4_helpers.R' and 'setMethod("print"' in first:
   target=('R/s4_helpers.R','_anonymous_.83');basis='same S4 print(sceptre_object) method'
  elif item['file']=='R/plotting_functions.R' and 'function(curr_grna_id)' in first:
   target=('R/plotting_functions.R','_anonymous_.154');basis='same grnas_to_plot lapply graph-construction site'
  else:
   target=('unmapped anonymous site',item['name']);basis='no inherited anonymous site established'
 matches=[x for x in old if (x['file'],x['name'])==target and int(x['length'])>50]
 inherited=bool(matches)
 y=matches[0] if inherited else {}
 lineage.append(dict(item,status='inherited_function_site' if inherited else 'NEW_LONG_FUNCTION',baseline_file=y.get('file',''),baseline_name=y.get('name',''),baseline_length=y.get('length',''),baseline_start_line=y.get('startLine',''),baseline_end_line=y.get('endLine',''),basis=basis))
if lineage:
 with (out/'function-lineage.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(lineage[0]));w.writeheader();w.writerows(lineage)
def long_lines(root):
 rows=[]
 files=sorted((root/'R').glob('*.R'))+sorted((root/'R').glob('*.r'))+sorted((root/'vignettes').glob('*.Rmd'))+[root/'NAMESPACE']
 for path in files:
  if not path.is_file():continue
  for i,text in enumerate(path.read_text().splitlines(),1):
   if len(text)>80:rows.append(dict(file=str(path.relative_to(root)),line=i,width=len(text),text=text))
 return rows
old_lines=long_lines(submitted);new_lines=long_lines(candidate);locations=collections.defaultdict(collections.deque)
for item in old_lines:locations[item['text']].append(item)
for item in new_lines:
 inherited=bool(locations[item['text']])
 item['status']='unchanged_submission_line' if inherited else 'NEW_OR_MODIFIED_LONG_LINE'
 baseline_item=locations[item['text']].popleft() if inherited else {}
 item['baseline_file']=baseline_item.get('file','')
 item['baseline_line']=baseline_item.get('line','')
 item['baseline_width']=baseline_item.get('width','')
if new_lines:
 with (out/'long-line-lineage.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(new_lines[0]));w.writeheader();w.writerows(new_lines)
result=dict(baseline_long_functions=sum(int(x['length'])>50 for x in old),candidate_long_functions=len(lineage),new_long_function_instances=sum(x['status']=='NEW_LONG_FUNCTION' for x in lineage),baseline_long_lines=len(old_lines),candidate_long_lines=len(new_lines),new_or_modified_long_line_instances=sum(x['status']=='NEW_OR_MODIFIED_LONG_LINE' for x in new_lines),scope='Exact submitted source and checker function-length algorithm; explicit reviewed semantic rename/anonymous-method lineage; long-line reuse requires exact source-text multiplicity. This is not an algorithmic parity assertion.')
result['passed']=result['new_long_function_instances']==0 and result['new_or_modified_long_line_instances']==0
(out/'instance-gate.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if not result['passed']:sys.exit(1)
