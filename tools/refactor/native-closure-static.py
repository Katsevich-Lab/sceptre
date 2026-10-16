#!/usr/bin/env python3
from pathlib import Path
import hashlib,importlib.metadata,json,re,subprocess
import lizard
if importlib.metadata.version("lizard") != "1.24.1":
 raise SystemExit("Source readback requires lizard 1.24.1.")
import argparse
parser=argparse.ArgumentParser(description="CL4 source-body parity and deletion readback; requires lizard 1.24.1; no R/compiled tests.")
parser.add_argument('root',type=Path)
parser.add_argument('output',type=Path)
parser.add_argument('--base',default='d3746f6f9aecd0dfb6bf954ce982dd81659b0026')
args=parser.parse_args()
root=args.root.resolve()
base=subprocess.check_output(['git','rev-parse',args.base],cwd=root,text=True).strip()
renames={'run_low_level_test_full_v4':'run_association_test','crt_index_sampler_fast':'sample_crt_indices','sample_combinations_v2':'sample_calibration_groups','sample_undercover_pairs_v2':'sample_calibration_pairs','increment_matrix':'increment_group_indices_in_place','legacy_result_record':'association_result_record','construct_negative_control_pairs_v2':'construct_negative_control_pairs'}
excluded={'RcppExports.cpp','compute_nb_size_param_functions.cpp','mixture_functs.cpp'}
def canonical(s):
 for a,b in renames.items():s=re.sub(r'\b'+a+r'\b',b,s)
 return s
comment=re.compile(r'("(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\')|(//[^\n]*|/\*.*?\*/)',re.S)
def tokens(s):
 s=comment.sub(lambda m:m.group(1) or ' ',s)
 return re.findall(r'[A-Za-z_]\w*|(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?|[^\s]',canonical(s))
def functions(path,source):
 lines=source.splitlines(True)
 fs={}
 for f in lizard.analyze_file.analyze_source_code(path,source).function_list:
  text=''.join(lines[f.start_line-1:f.end_line]);start=text.index('{');body=text[start:]
  name=canonical(f.name.split('::')[-1])
  fs[name]={'name':name,'start_line':f.start_line,'end_line':f.end_line,'signature_tokens':tokens(text[:start]),'body_tokens':tokens(body)}
 return fs
current={};baseline={}
paths=subprocess.check_output(['git','ls-tree','-r','--name-only',base,'--','src'],cwd=root,text=True).splitlines()
for relative in paths:
 p=Path(relative)
 if p.suffix not in ('.h','.cpp') or p.name in excluded:continue
 text=subprocess.check_output(['git','show',base+':'+relative],cwd=root,text=True)
 for name,f in functions(relative,text).items():baseline[(relative,name)]=f
for p in sorted((root/'src').iterdir()):
 if p.suffix not in ('.h','.cpp') or p.name in excluded:continue
 relative=p.relative_to(root).as_posix();text=p.read_text()
 for name,f in functions(relative,text).items():current[(relative,name)]=f
expected_removed={('src/check_functions.cpp',n) for n in ('synth_idx_list_to_matrix','synth_idx_list_to_r_list','print_synth_idx_list_row','test')}
expected_removed|={('src/low_level_full_test.cpp','compute_observed_full_statistic_v2'),('src/low_level_full_test.cpp','compute_null_full_statistics'),('src/inference_scores.cpp','estimate_log_fold_change_v2'),('src/generate_samples_functions.cpp','crt_index_sampler'),('src/sparse_matrix_functs.cpp','compute_colwise_max'),('src/sparse_matrix_functs.cpp','increment_vector')}
assert set(baseline)-set(current)==expected_removed,(set(baseline)-set(current))^expected_removed
assert not(set(current)-set(baseline))
for key,f in current.items():
 assert f['body_tokens']==baseline[key]['body_tokens'],key
 expected=baseline[key]['signature_tokens']
 if key==('src/inference_tail.cpp','check_for_outliers'):
  expected=expected.copy();expected.insert(expected.index('(')+1,'const')
 assert f['signature_tokens']==expected,('signature',key)
for relative in ('R/association_workers.R','R/s4_analysis_functs_2.R','R/neg_control_functions.R'):
 old=canonical(subprocess.check_output(['git','show',base+':'+relative],cwd=root,text=True))
 if relative=='R/neg_control_functions.R':
  old=old.replace('    grna_target_data_frame <- sceptre_object@grna_target_data_frame\n','').replace('    low_moi <- sceptre_object@low_moi\n','')
 assert old==(root/relative).read_text(),relative
frozen=subprocess.check_output(['git','diff','--name-only',base,'--','tests/testthat/fixtures/refactor','NAMESPACE','DESCRIPTION','R/storage.R','R/import_functs.R','src/compute_nb_size_param_functions.cpp','src/mixture_functs.cpp'],cwd=root,text=True)
assert not frozen,frozen
out={'base_commit':base,'scope':'Source-only structural check, not compiled/behavioral qualification or a complexity metric.',
 'retained_handwritten_implementations':len(current),'removed_implementations':len(expected_removed),
 'retained_body_tokens_equal_after_exact_symbol_crosswalk':True,
 'retained_signatures_equal_after_crosswalk_except_private_outlier_const':True,
 'R_callers_equal_after_crosswalk_and_two_unread_local_deletions':True,
 'frozen_and_excluded_source_untouched':True,
 'source_sha256':{p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((root/'src').iterdir()) if p.suffix in ('.cpp','.h') and p.name not in excluded},
 'functions':[dict(path=k[0],name=k[1],start_line=f['start_line'],end_line=f['end_line'],body_tokens_sha256=hashlib.sha256(json.dumps(f['body_tokens']).encode()).hexdigest(),signature_tokens_sha256=hashlib.sha256(json.dumps(f['signature_tokens']).encode()).hexdigest()) for k,f in sorted(current.items())],
 'removed':[dict(path=p,name=n) for p,n in sorted(expected_removed)]}
args.output.write_text(json.dumps(out,indent=2)+'\n')
print('PASS',len(current),'retained native bodies token-identical after approved renames;',len(expected_removed),'retirements; R callers canonical; protected sources unchanged.')
