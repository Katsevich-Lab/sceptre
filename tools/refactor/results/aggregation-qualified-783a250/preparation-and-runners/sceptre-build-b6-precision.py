from pathlib import Path
import subprocess,os,tarfile,sys,json
role=sys.argv[1]
repo=Path('/Users/ekatsevi/code/research/sceptre')
commit={'reference':'21f9ea098b69c444a884a49bc254de711968e3b5','candidate':'9631ef2'}[role]
pin=subprocess.check_output(['git','rev-parse',commit],cwd=repo,text=True).strip()
root=Path('/private/tmp/sceptre-b6-precision-9631ef2')/role
root.mkdir(parents=True,exist_ok=False)
source=root/'source'; library=root/'library'
source.mkdir();library.mkdir()
archive=root/'source.tar'
with archive.open('wb') as stream:
    subprocess.run(['git','archive','--format=tar',pin],cwd=repo,stdout=stream,check=True)
with tarfile.open(archive) as package:
    package.extractall(source,filter='data')
env=os.environ.copy(); env['R_LIBS']='/private/tmp/sceptre-rf0-ondisc-library'
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    env[k]='1'
log=root/'install.log'
with log.open('w') as stream:
    result=subprocess.run(['R','CMD','INSTALL','--preclean','--library='+str(library),str(source)],env=env,stdout=stream,stderr=subprocess.STDOUT)
if result.returncode: raise SystemExit('Install failed; inspect '+str(log))
subprocess.run(['python3',str(repo/'tools/refactor/resource-build-provenance.py'),str(repo),str(source),str(library),str(log),str(root/'build.json'),'--commit',pin],check=True)
print(json.dumps({'role':role,'pin':pin,'build':str(root/'build.json')}))
