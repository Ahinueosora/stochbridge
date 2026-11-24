"""
A6 step 1: reproduce Lamoureux 2023 Fig. 6 (anaerobic -> aerobic transition, project 'aat' from the K-12 modulome)
with our re-implementation (src/dima.py) and compare against the table saved in the authors' notebook.

Writes outputs/a6/aat_dima_reproduction.csv and outputs/a6/aat_check.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from dima import P1K_DIR, KnowledgeBase, log_tpm_from_counts  # noqa: E402

ROOT = Path(__file__).parent.parent
OUT = ROOT / "outputs" / "a6"

# DiMA table saved in precise1k/notebooks/precise1k/workflow__add_new_data.ipynb (fdr 0.05, threshold 9)
PAPER_AAT = {"ArcA": 26.53, "Fur-2": 17.51, "Fur-1": 13.87, "cyoB/kdpE/qseB KO": 13.35, "sdh KO": 9.95, "IS Elements-1": 9.18,
             "FDH-N/NRA": -9.08, "Cytochrome c": -10.01, "Microaerobic": -10.01, "Fnr-2": -11.48, "Fnr-1": -14.72}


def gene_lengths(kb):
    return (kb.gene_table.end - kb.gene_table.start).abs() + 1


def load_k12(kb, projects):
    md = pd.read_csv(P1K_DIR / "k12_metadata_qc.csv", index_col=0, low_memory=False)
    md = md[md.project.isin(projects)]
    first = pd.read_csv(P1K_DIR / "k12_counts.csv", nrows=0).columns[0]
    keep = set(md.index) | {first}
    counts = pd.read_csv(P1K_DIR / "k12_counts.csv", usecols=lambda c: c in keep).set_index(first)
    md = md.loc[md.index.intersection(counts.columns)]
    return log_tpm_from_counts(counts[md.index], gene_lengths(kb)), md


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    kb = KnowledgeBase.load()
    lt, md = load_k12(kb, ["aat"])
    kb.add_samples(lt, md)
    res = kb.dima("aat:00m_00s", "aat:10m_00s", threshold=9, fdr=0.05)
    res.to_csv(OUT / "aat_dima_reproduction.csv")
    sig = res[res.significant]
    paper = pd.Series(PAPER_AAT)
    both = paper.index.intersection(res.index)
    chk = {"n_samples": int(len(md)), "ours_significant": sig.difference.round(2).to_dict(),
           "paper_significant": PAPER_AAT,
           "overlap": sorted(set(sig.index) & set(paper.index)), "only_ours": sorted(set(sig.index) - set(paper.index)),
           "only_paper": sorted(set(paper.index) - set(sig.index)),
           "difference_r_vs_paper": float(pd.concat([res.loc[both, "difference"], paper[both]], axis=1).corr().iloc[0, 1]),
           "max_abs_diff_vs_paper": float((res.loc[both, "difference"] - paper[both]).abs().max())}
    (OUT / "aat_check.json").write_text(json.dumps(chk, indent=2))
    return chk


if __name__ == "__main__":
    print(json.dumps(main(), indent=1))
