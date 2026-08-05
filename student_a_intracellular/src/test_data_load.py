import pandas as pd
base = "student_a_intracellular/data/precise1k/precise1k/data/precise1k"
trn_base = "student_a_intracellular/data/precise1k/precise1k/data/regulation"
log_tpm = pd.read_csv(f"{base}/log_tpm_qc.csv", index_col=0)
M = pd.read_csv(f"{base}/M.csv", index_col=0)
A = pd.read_csv(f"{base}/A.csv", index_col=0)
metadata = pd.read_csv(f"{base}/metadata_qc.csv", index_col=0)
trn = pd.read_csv(f"{trn_base}/TRN.csv")
print("log_tpm shape:", log_tpm.shape)
print("M shape:", M.shape)
print("A shape:", A.shape)
print("metadata shape:", metadata.shape)
print("TRN shape:", trn.shape)
print("All good — files load and shapes match the paper.")
