"""
Download the input data for Student A's pipeline (run once, from student_a_intracellular/).

    python get_data.py

What it does
  1. Fetches the 10 PRECISE-1K / K-12 modulome files the pipeline needs from github.com/SBRG/precise1k,
     pinned to commit 4829b83 (29 Aug 2023), into data/precise1k/precise1k/ (about 240 MB).
     Only those files are downloaded (git partial clone + sparse checkout), not the whole repository.
  2. Checks every file against its SHA-256 checksum, so you know you have exactly the data the results were made from.
  3. Makes the short names the code reads (data/precise1k/p1k_*.csv, k12_*.csv): symbolic links, or copies where
     links are not allowed (some Windows set-ups).

The PETRI-seq single-cell data (GEO GSE141018, sample GSM4489500, 3.9 MB) is small and is already in the
repository under data/scrnaseq/GSE141018/, so nothing needs downloading for A5.

Needs: git (any recent version).
"""
import hashlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO = "https://github.com/SBRG/precise1k"
COMMIT = "4829b83ace51ee7980de6950077f1cf8f3746f88"
ROOT = Path(__file__).resolve().parent
DEST = ROOT / "data" / "precise1k" / "precise1k"

FILES = {  # path in SBRG/precise1k -> SHA-256
    "data/precise1k/M.csv": "6ad63e36378b9da8e7c5b501cca81a5db7f65872684f4a870428e37380c06afa",
    "data/precise1k/A.csv": "a91eea1b57c14c94aea379d40496ada6e0af32f9e708cbc51e0329ea570608f7",
    "data/precise1k/metadata_qc.csv": "a1358525a0431ac68b44f85f56169b7e8af4efb9186be6b43ba2cacb7b56e1cf",
    "data/precise1k/counts.csv": "f6307411d59e640d46448353e558c5761bf2b95e20f3d098bbbc10af93e95e10",
    "data/precise1k/log_tpm_qc.csv": "8bdb286aecd305a38323ad5a5c18eb264702739824dfcb708c5c5b595abb84d6",
    "data/precise1k/precise1k.json.gz": "2f5f0e41867fc2121cefdde8aeb4efd545466d35d6ad42b604d402db62d67de9",
    "data/k12_modulome/counts.csv": "e60ac8370095e180badc89179a1cfd448c748c2ebf9dde0a31d2dc381318ec20",
    "data/k12_modulome/metadata_qc.csv": "3904fb6db7c30396899dd78b045e0b18358ae0f83dbb33e6984024532e84d53e",
    "data/regulation/TRN.csv": "d30578836279465414522061ccc9b4eae504b7df36725380491c39caec9b4089",
    "data/annotation/gene_info.csv": "355de1157a4837180cb0c3d625a7a4f26903ba0d50f2bf3ca769949e90dbcb90",
}

SHORT_NAMES = {  # name the code reads (in data/precise1k/) -> file in the clone
    "p1k_precise1k.json.gz": "data/precise1k/precise1k.json.gz",
    "p1k_log_tpm_qc.csv": "data/precise1k/log_tpm_qc.csv",
    "p1k_metadata_qc.csv": "data/precise1k/metadata_qc.csv",
    "p1k_counts.csv": "data/precise1k/counts.csv",
    "k12_counts.csv": "data/k12_modulome/counts.csv",
    "k12_metadata_qc.csv": "data/k12_modulome/metadata_qc.csv",
}


def git(*args, cwd=None):
    subprocess.run(["git", "-c", "gc.auto=0", "-c", "advice.detachedHead=false", *args], cwd=cwd, check=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fetch():
    if not (DEST / ".git").exists():
        DEST.parent.mkdir(parents=True, exist_ok=True)
        print(f"Cloning {REPO} (only the needed files) into {DEST.relative_to(ROOT)} ...")
        git("clone", "--filter=blob:none", "--no-checkout", REPO, str(DEST))
    git("sparse-checkout", "set", "--no-cone", *FILES, cwd=DEST)
    git("checkout", "--quiet", COMMIT, cwd=DEST)


def check():
    bad = []
    for rel, want in FILES.items():
        p = DEST / rel
        ok = p.exists() and sha256(p) == want
        print(f"  {'OK ' if ok else 'BAD'} {rel}")
        if not ok:
            bad.append(rel)
    if bad:
        sys.exit(f"{len(bad)} file(s) missing or different from the pinned version: {bad}")


def link():
    base = ROOT / "data" / "precise1k"
    for short, rel in SHORT_NAMES.items():
        dst, src = base / short, DEST / rel
        if dst.is_symlink() or dst.exists():
            dst.unlink()
        try:
            os.symlink(os.path.relpath(src, base), dst)
        except OSError:
            shutil.copyfile(src, dst)


if __name__ == "__main__":
    fetch()
    print("Checking files against the pinned checksums:")
    check()
    link()
    petri = ROOT / "data" / "scrnaseq" / "GSE141018"
    n = len(list(petri.glob("GSM4489500_*.txt.gz"))) if petri.exists() else 0
    print("PETRI-seq data (in the repository):", "OK" if n else "MISSING - re-clone the repository")
    print("Done. Data is ready in data/.")
