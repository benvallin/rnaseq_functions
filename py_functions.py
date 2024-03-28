# %% Set up ----

# Import required libraries
import numpy as np
import pandas as pd
import scanpy as sc

# %% find_top_expressing_cells() ----

def find_top_expressing_cells(adata, gene_name: str, top_pct: int):
  
  if gene_name not in list(adata.var["gene_name"]):
    return print(f"{gene_name} not in adata.var[\"gene_name\"]")
  
  top_quantile = 1-(top_pct/100)

  gene_counts = adata[:, adata.var["gene_name"]==gene_name].X[:, 0].copy()

  top_quantile_value = np.quantile(a=gene_counts, q=top_quantile, axis=0)
  
  top_pct_expressing = np.where(gene_counts >= top_quantile_value, True, False)
  
  top_pct_expressing = pd.Series(data=top_pct_expressing, index=adata.obs_names, dtype="category")

  if len(set(top_pct_expressing)) > 1:
    top_pct_expressing = top_pct_expressing.cat.reorder_categories(new_categories=[True, False])
    adata.obs["top"+str(top_pct)+"pct_"+gene_name+"_expressing"] = top_pct_expressing
    adata.uns["top"+str(top_pct)+"pct_"+gene_name+"_expressing_colors"] = ["#4361ee", "#e7ece5"]
  else:
    adata.obs["top"+str(top_pct)+"pct_"+gene_name+"_expressing"] = top_pct_expressing
    adata.uns["top"+str(top_pct)+"pct_"+gene_name+"_expressing_colors"] = ["#907ad6"]

# %% find_top_expressing_cells() ----
    
def find_least_expressing_cells(adata, gene_name: str, least_pct: int):
  
  if gene_name not in list(adata.var["gene_name"]):
    return print(f"{gene_name} not in adata.var[\"gene_name\"]")
  
  least_quantile = 1-(least_pct/100)

  gene_counts = adata[:, adata.var["gene_name"]==gene_name].X[:, 0].copy()

  least_quantile_value = np.quantile(a=gene_counts, q=least_quantile, axis=0)
  
  least_pct_expressing = np.where(gene_counts <= least_quantile_value, True, False)
  
  least_pct_expressing = pd.Series(data=least_pct_expressing, index=adata.obs_names, dtype="category")

  if len(set(least_pct_expressing)) > 1:
    least_pct_expressing = least_pct_expressing.cat.reorder_categories(new_categories=[True, False])
    adata.obs["least"+str(least_pct)+"pct_"+gene_name+"_expressing"] = least_pct_expressing
    adata.uns["least"+str(least_pct)+"pct_"+gene_name+"_expressing_colors"] = ["#f72585", "#e7ece5"]
  else:
    adata.obs["least"+str(least_pct)+"pct_"+gene_name+"_expressing"] = least_pct_expressing
    adata.uns["least"+str(least_pct)+"pct_"+gene_name+"_expressing_colors"] = ["#907ad6"]