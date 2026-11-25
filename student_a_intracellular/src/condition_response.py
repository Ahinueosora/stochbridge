"""
Extracts a regulon's iModulon activity as a function of a chosen environmental
or genetic variable, and flags samples matching a strain-description keyword
(e.g. CRP genetic perturbations).
"""

import pandas as pd


def categorize_by_strain_keyword(sample_table, keyword, column="Strain Description"):
    """
    Flags samples whose strain description contains a keyword (case-insensitive).

    Parameters
    ----------
    sample_table : pd.DataFrame
    keyword : str
        e.g. "crp"
    column : str
        Column in sample_table to search.

    Returns
    -------
    pd.Series (bool), indexed like sample_table
    """
    return sample_table[column].astype(str).str.contains(keyword, case=False, na=False)


def categorize_crp_strains(strain_description):
    """
    Categorizes a single CRP-related strain description string into one of:
    'CRP full deletion', 'CRP ΔAr1 (Class I)', 'CRP ΔAr2 (Class II)',
    'CRP ΔAr1+Ar2', or 'Other/WT'.

    Parameters
    ----------
    strain_description : str

    Returns
    -------
    str
    """
    s = str(strain_description).lower().replace(" ", "")
    if "delar1delar2" in s:
        return "CRP ΔAr1+Ar2"
    elif "delar2" in s:
        return "CRP ΔAr2 (Class II)"
    elif "delar1" in s:
        return "CRP ΔAr1 (Class I)"
    elif "del_crp" in s or "delcrp" in s:
        return "CRP full deletion"
    else:
        return "Other/WT"


def get_condition_response(ica_data, imodulon_x, imodulon_y, category_fn=None,
                             category_column="Strain Description"):
    """
    Builds a DataFrame of two iModulons' activities per sample, optionally
    with a category label per sample (e.g. genetic perturbation type).

    Parameters
    ----------
    ica_data : imodulon_loader.IcaData
    imodulon_x, imodulon_y : str
        Names of the two iModulons to compare, e.g. "Crp-1", "Crp-2".
    category_fn : callable, optional
        Function applied to sample_table[category_column] to produce a
        category label per sample.
    category_column : str
        Column in sample_table to pass to category_fn.

    Returns
    -------
    pd.DataFrame
        Columns: [imodulon_x, imodulon_y, 'category' (if category_fn given)]
    """
    df = pd.DataFrame({
        imodulon_x: ica_data.A.loc[imodulon_x],
        imodulon_y: ica_data.A.loc[imodulon_y],
    })

    if category_fn is not None:
        df["category"] = ica_data.sample_table[category_column].apply(category_fn)

    return df
