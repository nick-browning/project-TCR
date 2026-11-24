import pandas as pd

REF_BARCODE_LIST = "/data/projects/punim1612/CSI/Nick-project/Slurm_scripts/CITE_barcode_list.txt"
BARCODE_LIST=[]

with open(REF_BARCODE_LIST, "r") as f:
    for BARCODE in f:
        BARCODE_LIST.append(BARCODE.strip())

reference = pd.read_csv("totalseq-C_antibodies.csv")
barcodes = reference[reference["id"].isin(BARCODE_LIST)]

barcodes.to_csv("cellranger_multi_feature_ref.csv", index=False)