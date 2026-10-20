#!/usr/bin/env python3
"""Build a frozen Git revision, or qualify an exact supplied source archive.
No dependency installation. Every install targets output/library. Output must be
new so old pass receipts cannot survive a failed rerun. Invoke under local lease.
"""
import argparse, hashlib, json, os, pathlib, shutil, subprocess, sys, tarfile, time

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=pathlib.Path, default=pathlib.Path.cwd())
    p.add_argument("--revision", default="HEAD")
    p.add_argument("--archive", type=pathlib.Path)
    p.add_argument("--output", type=pathlib.Path, required=True)
    p.add_argument("--ondisc", choices=["present", "absent", "either"], default="either")
    p.add_argument("--check", action="store_true")
    p.add_argument("--as-cran", action="store_true",
                   help="Use hosted --as-cran checks including vignette rebuilds")
    p.add_argument("--r", default="R")
    p.add_argument("--rscript", default="Rscript")
    a = p.parse_args()
    if a.as_cran and not a.check:
        p.error("--as-cran requires --check")
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ)
    env.update({x: "1" for x in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"]})
    # Existing dependency libraries are read-only fallbacks, never install targets.
    dependencies = subprocess.check_output([a.rscript, "--vanilla", "-e",
        'cat(.libPaths(), sep="\\n")'], text=True, env=env).splitlines()
    env["R_LIBS"] = os.pathsep.join(dependencies)
    env.update(MAKEFLAGS="-j2", R_PROFILE_USER=os.devnull, R_ENVIRON_USER=os.devnull,
               R_LIBS_USER=str(out / "library"))
    (out / "library").mkdir()
    receipt = dict(status="running", started_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                   commands=[], ondisc_policy=a.ondisc)
    def save():
        (out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    def run(cmd, log, cwd=out):
        entry = dict(argv=list(map(str, cmd)), cwd=str(cwd), log=log)
        receipt["commands"].append(entry)
        save()
        with (out / log).open("w") as f:
            proc = subprocess.run(entry["argv"], cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT)
        entry["exit_code"] = proc.returncode
        save()
        if proc.returncode:
            raise RuntimeError(f"{log}: exit {proc.returncode}")
    try:
        if a.archive:
            archive = out / a.archive.name
            shutil.copyfile(a.archive.resolve(), archive)
            receipt["input_archive"] = str(a.archive.resolve())
        else:
            source = a.source.resolve()
            revision = subprocess.check_output(["git", "-C", str(source), "rev-parse",
                a.revision + "^{commit}"], text=True).strip()
            receipt["source_commit"] = revision
            receipt["source_tree"] = subprocess.check_output(["git", "-C", str(source),
                "rev-parse", revision + "^{tree}"], text=True).strip()
            run(["git", "-C", source, "archive", "--format=tar", "--prefix=sceptre/",
                "-o", out / "git-source.tar", revision], "git-archive.log")
            with tarfile.open(out / "git-source.tar") as t:
                # git archive supplies repository paths, not an untrusted download.
                t.extractall(out / "build-source")
            run([a.r, "CMD", "build", out / "build-source/sceptre", "--no-manual"], "build.log")
            archives = list(out.glob("sceptre_*.tar.gz"))
            if len(archives) != 1:
                raise RuntimeError("Expected exactly one built archive")
            archive = archives[0]
        receipt["archive"] = archive.name
        receipt["archive_sha256"] = hashlib.sha256(archive.read_bytes()).hexdigest()
        save()
        run([a.r, "CMD", "INSTALL", "--preclean", "--install-tests",
            "--library=" + str(out / "library"), archive], "install.log")
        run([a.rscript, "--vanilla", pathlib.Path(__file__).with_name("installed-suite.R"),
            out / "library", out, a.ondisc], "tests.log")
        if a.check:
            check_args = ["--no-manual", "--as-cran"] if a.as_cran else [
                "--no-manual", "--no-build-vignettes"]
            run([a.r, "CMD", "check", archive] + check_args, "check.log")
            check = (out / "sceptre.Rcheck/00check.log").read_text()
            receipt["R_CMD_check_status"] = [x for x in check.splitlines() if x.startswith("Status:")]
            if not receipt["R_CMD_check_status"]:
                raise RuntimeError("Missing terminal R CMD check status")
            if any(x in " ".join(receipt["R_CMD_check_status"]) for x in ["ERROR", "WARNING"]):
                raise RuntimeError("R CMD check errors/warnings; inspect retained diagnostics")
        receipt["status"] = "passed"
    except Exception as e:
        receipt.update(status="failed", error=str(e))
    finally:
        receipt["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        save()
    print(json.dumps(receipt, indent=2))
    return 0 if receipt["status"] == "passed" else 1

if __name__ == "__main__":
    sys.exit(main())
