import sys, scanpy as sc

selected_dataset = [
    'Kidney_BenPublished', 
    'Kidney_BenUnpublished',
    'Kidney_Krebs',
    'Kidney_Liao',
    'Kidney_Malone',
    'Kidney_Muto', 
    'Kidney_Raji', 
    'Kidney_Wilson',
    'Kidney_Krishna',
    'Kidney_Wu'
]

# https://docs.google.com/spreadsheets/d/1yS3yNlVnIlsGclgao0R9EeTPQTKF8Nx4wPupjzM3mrg/edit?gid=0#gid=0

cell_types = [
    "SULF1+ EC-AEA", "EC-AEA", "Macro", "Macro", "cycMacro", "CCD-IC-A",
    "OMCD-IC-A", "CNT-PC", "CNT-PC", "DTL", "DTL", "DTL", "EC-DVR",
    "SLC6A6+ EC-AEA", "EC-GC", "EC-PTC", "EC-PTC", "EC-AVR", "IC-B",
    "VSMC", "Pericyte", "MAST", "MC", "MD", "cycNK/T", "cycNK/T", "DC1",
    "DC2", "IC-B|CNT doub", "POD", "POD", "PEC", "PT-S1/2", "PT-S1/2_nuc",
    "PT-S3", "dPT", "dPT", "TAL", "TAL", "TAL", "ATL", "Treg", "CD4 T",
    "CD4 T", "CD4 T", "Th17", "CD8 T", "CD8 T", "CD8 T", "ILC3", "NKT",
    "Ciliated", "DCT", "EC-LYM", "OMCD-PC", "IMCD-PC", "MFAP5+aFIB",
    "aFIB", "aFIB", "aFIB", "aPT", "aPT", "CCD-PC", "CCD-PC", "DCT",
    "Neuron", "RBC", "B", "B", "PL", "TAL", "VWF+ EC-AVR", "VWF+ EC-AVR",
    "cycEC", "cycPT", "CD16+ NK", "CD56bright NK", "pDC", "cMono", "ncMono",
    "cycPapE", "PapE", "γδT", "dC-IC-A"
]



if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python labeled_subset.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    print(f"Reading {input_file}")
    adata = sc.read_h5ad(input_file)
    adata_labeled = adata[adata.obs["harmonized_celltype"] != 'nan']
    adata_labeled = adata_labeled[[_batch in selected_dataset for _batch in adata_labeled.obs["Project_ID"]]]
    adata_harmonized = adata_labeled[[celltype in cell_types for celltype in adata_labeled.obs["harmonized_celltype"]]]
    
    if "scPoli" in adata_harmonized.obsm.keys():
        del adata_harmonized.obsm["scPoli"]
    
    print(f"Saving to {output_file}")
    adata_harmonized.write_h5ad(output_file, compression="gzip")
    