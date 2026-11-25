# Curated inputs (small, tracked in git)

| File | What it is | Source | Used by |
|---|---|---|---|
| `regulondb_crp_edges.json` | RegulonDB evidence for each CRP → target edge in the circuit (effect, confidence level, evidence type, PubMed IDs, causal grade 1–4) | RegulonDB release 14.5.1, TF-RISet download from regulondb.ccg.unam.mx, accessed 24 Nov 2026 | `src/build_a3_v3.py` (A3) |
| `lac_network_gene_check.xlsx` | The combined CRP/lac network table: RegulonDB edges, the 21 A2 causal edges, Kuhlman 2007 and Ozbudak 2004 edges, with evidence | Assembled by hand from the sources above (sheet "Combined network") | `src/build_network.py` (A3) |

RegulonDB is updated over time, so these files are kept here rather than downloaded, to keep the results reproducible.
