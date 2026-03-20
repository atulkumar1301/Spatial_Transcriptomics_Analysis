import scanpy as sc
import pandas as pd

# Load your .h5ad single-cell data
adata = sc.read_h5ad("/Users/akumar/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Projects/Minna_Spatial_Transcriptomics_Data/8f884f96-193f-429d-a6c5-13d68c192ba5.h5ad")

# Load the mouse TSV mapping file

mapping_df = pd.read_csv("/Users/akumar/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Projects/Minna_Spatial_Transcriptomics_Data/Gene_Mapping_Ensembel_to_Gene_Name.txt", sep='\t')
mapping_dict = dict(zip(mapping_df['ensembl_gene_id'], mapping_df['mgi_symbol']))

# Backup original IDs and map to new symbols
# This assumes the current index (adata.var_names) contains the Ensembl IDs
adata.var['ensembl_gene_id'] = adata.var_names
adata.var['mgi_symbol'] = adata.var_names.map(mapping_dict)

# Handle missing mappings
# If an ID isn't in your TSV, keep the original Ensembl ID to avoid 'NaN' values
adata.var['mgi_symbol'] = adata.var['mgi_symbol'].fillna(adata.var['ensembl_gene_id'])

# Update the index to Gene Symbols
adata.var_names = adata.var['mgi_symbol'].astype(str)

# Resolve duplicates
# Multiple Ensembl IDs sometimes map to one symbol; Scanpy requires unique names
adata.var_names_make_unique()

# Save the updated object
adata.write('/Users/akumar/Library/CloudStorage/OneDrive-UniversityofEasternFinland/Projects/Minna_Spatial_Transcriptomics_Data/Updated_mouse_Aortas_data.h5ad')
