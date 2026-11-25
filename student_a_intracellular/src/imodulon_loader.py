"""
Loads PRECISE-1K M, A, sample_table, and TRN matrices into a small IcaData container.
Handles known data-quirks in the SBRG/precise1k repo:
- Must use *_qc.csv files (log_tpm_qc, metadata_qc), not the raw versions,
  to match the 4257 gene x 1035 sample dimensions reported in the paper.
- 'Superoxide' in M.csv and 'SoxS' in A.csv refer to the same iModulon
  (position 22 in both) but are named differently across files; we rename
  M's column to match A's for consistency.

IcaData below holds the same matrices as pymodulon.core.IcaData and checks the same thing
(M and A must name the same iModulons). It replaces pymodulon here because the published
pymodulon 0.2.1 does not install or import with current pandas; A1 only needs the matrices.
Data: run `python get_data.py` once (puts the SBRG/precise1k files in data/precise1k/precise1k/).
"""

from pathlib import Path

import pandas as pd

_DATA = Path(__file__).resolve().parents[1] / "data" / "precise1k" / "precise1k" / "data"
DEFAULT_BASE = str(_DATA / "precise1k")
DEFAULT_TRN_BASE = str(_DATA / "regulation")


class IcaData:
    """Container for an ICA decomposition: M (genes x iModulons), A (iModulons x samples), metadata, TRN."""

    def __init__(self, M, A, gene_table=None, sample_table=None, trn=None):
        if list(M.columns) != list(A.index):
            raise ValueError("M and A matrices have different iModulon names")
        self.M, self.A = M, A
        self.gene_table, self.sample_table, self.trn = gene_table, sample_table, trn

    @property
    def imodulon_names(self):
        return list(self.A.index)


def load_precise1k_matrices(base=DEFAULT_BASE, trn_base=DEFAULT_TRN_BASE):
    """
    Loads the raw M, A, sample_table, and TRN matrices from the precise1k repo.

    Returns
    -------
    M, A, sample_table, trn : pd.DataFrame
    """
    M = pd.read_csv(f"{base}/M.csv", index_col=0)
    A = pd.read_csv(f"{base}/A.csv", index_col=0)
    sample_table = pd.read_csv(f"{base}/metadata_qc.csv", index_col=0)
    trn = pd.read_csv(f"{trn_base}/TRN.csv")

    # Known naming inconsistency: 'Superoxide' (M) vs 'SoxS' (A), same iModulon
    if "Superoxide" in M.columns and "SoxS" in A.index:
        M = M.rename(columns={"Superoxide": "SoxS"})

    return M, A, sample_table, trn


def build_ica_data(base=DEFAULT_BASE, trn_base=DEFAULT_TRN_BASE):
    """
    Loads matrices and builds a pymodulon IcaData object.

    Returns
    -------
    ica_data : IcaData
    """
    M, A, sample_table, trn = load_precise1k_matrices(base, trn_base)
    ica_data = IcaData(M, A, gene_table=None, sample_table=sample_table, trn=trn)
    return ica_data


def get_activity_by_condition(ica_data, imodulon_names):
    """
    Extracts activity values for one or more iModulons, joined with sample metadata.

    Parameters
    ----------
    ica_data : IcaData
    imodulon_names : list of str
        e.g. ["Crp-1", "Crp-2"]

    Returns
    -------
    pd.DataFrame
        Index: sample_id. Columns: one per iModulon, plus all sample_table columns.
    """
    activity_df = ica_data.A.loc[imodulon_names].T
    return activity_df.join(ica_data.sample_table)
