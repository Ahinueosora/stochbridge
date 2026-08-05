import pandas as pd
base = "student_a_intracellular/data/precise1k/precise1k/data/precise1k"
meta_qc = pd.read_csv(f"{base}/metadata_qc.csv", index_col=0)
print("metadata_qc shape:", meta_qc.shape)
