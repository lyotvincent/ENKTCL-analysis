from Bio import SeqIO
import pandas as pd
import scanpy as sc
import anndata as ad
# 聚类/重聚类代码
def clu(adata, key_added="majorType-fix", 
        n_neighbors=50, n_pcs=30, 
        rep='X_pca_harmony',  # X_pca
        do_har=False, # if do harmony
        max_iter=20, # harmony max iteration
        do_scrublet=False, # remove doublet cells 
        har_key='batch',   # harmony group obs
        resolution=1):  # cluster resolution
    # Computing the neighborhood graph
    if do_scrublet:
        n0 = adata.shape[0]
        print("{0} Cell number: {1}".format(key_added, n0))
        sc.external.pp.scrublet(adata)
        adata = adata[adata.obs['predicted_doublet']==False,:].copy()
        print("{0} Cells retained after scrublet, {1} cells reomved.".format(adata.shape[0], n0-adata.shape[0]))
    else:
        print("Ignoring processing doublet cells...")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.pca(adata, svd_solver='arpack', use_highly_variable=True)
    if do_har and len(adata.obs[har_key].cat.categories) > 1:
        sc.external.pp.harmony_integrate(adata, key=har_key,max_iter_harmony=max_iter)
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=rep)
    else:
        print("Evaluating neighbors only...")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=rep)
    # Run UMAP
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution, key_added=key_added)
    sc.pl.umap(adata, color=key_added, legend_fontoutline=True, palette=sc.pl.palettes.default_20, legend_loc="on data")
    return adata

def anno(adata:ad.AnnData, annoDict:dict, obsKey='cnv_status', obsVal='cnv_leiden', default="Unknown"):
  if default is not None:
    adata.obs[obsKey] = default
  for key in annoDict.keys():
    adata.obs.loc[adata.obs[obsVal].isin(annoDict[key]), obsKey] = key
  return adata


# memento差异分析
import memento
def DE_memento_1d(adata: ad.AnnData, ctleft: list, ctright: list, groupby='annotation'):
  """
  使用 Memento 包分析单细胞测序数据中的动态变化。

  Args:
      adata (ad.AnnData): 单细胞测序数据的 AnnData 对象。
      ctleft (list): 左条件列表，表示刺激前的细胞类型或状态。
      ctright (list): 右条件列表，表示刺激后的细胞类型或状态。
      groupby (str, optional): 用于分组的注释列名。默认为 'annotation'。

  Returns:
      pd.DataFrame: Memento 分析的结果，包括差异表达分析、HT 1D moments 分析等。

  """
  _adata = adata[adata.obs[groupby].isin(ctleft + ctright)].copy()
  _adata.obs['stim'] = _adata.obs[groupby].isin(ctright)
  _adata.obs['stim'] = _adata.obs['stim'].astype('int')
  _adata.obs['capture_rate'] = 0.07
  memento.setup_memento(_adata, q_column='capture_rate')
  memento.create_groups(_adata, label_columns=['stim'])
  memento.compute_1d_moments(_adata, min_perc_group=.7)
  meta_df = memento.get_groups(_adata)
  meta_df['intercept'] = 1
  convariate = meta_df[['intercept']]
  treatment = (meta_df[['stim']] == 1).astype(float)
  
  memento.ht_1d_moments(_adata,
  covariate=convariate,
  treatment=treatment,
  num_cpus=16,
  resample_rep=False,
  )
  memento_results = memento.get_1d_ht_result(_adata)
  return(memento_results)

# CytoTRACE2
import datatable as dt
from cytotrace2_py.cytotrace2_py import cytotrace2
def prep_cyto(adata:ad.AnnData, anno:str, rawpath:str, annopath:str):
  """
  将AnnData对象转换为CytoTrace2格式并写入文件。

  Args:
      adata (ad.AnnData): AnnData对象，包含待转换的数据。
      anno (str): 注释列的名称，用于生成包含细胞类型信息的文件。
      rawpath (str): 存储CytoTrace2格式数据的文件路径。
      annopath (str): 存储包含细胞类型信息的文件的路径。

  Returns:
      None

  """
  print("Converting AnnData to CytoTrace2 format...")
  cyto = adata.raw.X.astype('int')  # cytotrace2要求数据为原始数据，可以自行替换，默认是稀疏矩阵格式
  adata.obs_names_make_unique()
  adata.var_names_make_unique()
  cyto = pd.DataFrame(cyto.todense(), index=list(adata.obs_names), columns=list(adata.var_names)).T
  cyto = cyto.reset_index()
  cyto = dt.Frame(cyto)
  print("Writing CytoTrace2 files...")
  cyto.to_csv(rawpath, sep='\t', compression='gzip')
  cyto_anno = pd.DataFrame(list(adata.obs_names), columns=['cell IDs'])
  cyto_anno['phenotype'] = list(adata.obs[anno])
  cyto_anno.to_csv(annopath, sep='\t', index=False, header=True)
  print("Writing CytoTrace2 files finished...")

# KEGG、GO富集分析
import scanpy as sc
import gseapy as gp
from gseapy import Msigdb
def prep_enrichr(adata, groupby):
  # de analysis
  sc.tl.rank_genes_groups(adata, groupby=groupby)
  sc.tl.dendrogram(adata, groupby=groupby)
  # get deg result
  result = adata.uns['rank_genes_groups']
  groups = result['names'].dtype.names
  degs = pd.DataFrame(
      {group + '_' + key: result[key][group]
      for group in groups for key in ['names','scores', 'pvals','pvals_adj','logfoldchanges']})
  return degs

def run_enrichr_genelist(genelist):
  kegg = gp.enrichr(genelist,
                  gene_sets='KEGG_2021_Human',
                  outdir=None)
  cc = gp.enrichr(genelist,
                      gene_sets='GO_Cellular_Component_2023',
                      outdir=None)
  hm = gp.enrichr(genelist,
                      gene_sets=['MSigDB_Hallmark_2020'],
                      outdir=None)
  mf = gp.enrichr(genelist,
                      gene_sets='GO_Molecular_Function_2023',
                      outdir=None)
  bp = gp.enrichr(genelist,
                      gene_sets=['GO_Biological_Process_2023'],
                      outdir=None)
  return kegg, cc, hm, mf, bp

def run_enrichr_scanpy(degs, group,direction='up', top_term=15, all_term=False, plot=True):
  # subset up or down regulated genes
  degs_sig = degs[degs[f'{group}_pvals_adj'] < 0.05]
  degs_up = degs_sig[degs_sig[f'{group}_logfoldchanges'] > 0]
  degs_dw = degs_sig[degs_sig[f'{group}_logfoldchanges'] < 0]
  degs_used = degs_up[f'{group}_names'] if direction == 'up' else degs_dw[f'{group}_names']
  kegg, cc, hm, mf, bp = run_enrichr_genelist(degs_used)
  go_res = pd.concat([kegg.res2d.head(top_term),
                    hm.res2d.head(top_term), 
                    cc.res2d.head(top_term), 
                    bp.res2d.head(top_term), 
                    mf.res2d.head(top_term)])
  if plot:
    ax = gp.barplot(go_res, figsize=(6,20),
                group ='Gene_set',
                top_term=top_term,
                color = ['#8b9cc4', '#f08961', '#62bb9f', '#e71f19', '#3d5eaa'],
                title ="The Most Enriched Terms")
  else:
    ax = None
  if all_term:
    go_res = pd.concat([kegg.res2d, hm.res2d, cc.res2d, bp.res2d, mf.res2d])
  return go_res, ax
  
def run_enrichr_memento(memento_results,direction='up', top_term=15):
  memento_sig = memento_results[memento_results['de_pval'] < 0.05]
  memento_up = memento_sig[memento_sig['de_coef'] > 0]
  memento_dw = memento_sig[memento_sig['de_coef'] < 0]
  degs_used = memento_up['gene'] if direction == 'up' else memento_dw['gene']
  kegg, cc, hm, mf, bp = run_enrichr_genelist(degs_used)
  go_res = pd.concat([kegg.res2d.head(15),
                      hm.res2d.head(15), 
                      cc.res2d.head(15), 
                      bp.res2d.head(15), 
                      mf.res2d.head(15)])
  ax = gp.barplot(go_res, figsize=(6,20),
              group ='Gene_set',
              top_term=top_term,
              color = ['#8b9cc4', '#f08961', '#62bb9f', '#e71f19', '#3d5eaa'],
              title ="The Most Enriched Terms")
  return go_res, ax

# 细胞比例出图
import numpy as np
import matplotlib
from matplotlib.ticker import FuncFormatter
def countPlot(age_df, barlabels, ax, 
              colors=sc.pl.palettes.default_20, 
              xlabel='Patient ID', bar_label='total'):
    all_df = age_df.sum(axis=1)
    age_normdf = pd.DataFrame([age_df.loc[i,:] for i in barlabels], index=barlabels)
    age_norm1df = pd.DataFrame([age_df.loc[i,:]/all_df[i] for i in barlabels], index=barlabels)
    age_cumdf = pd.DataFrame([np.cumsum(age_normdf.loc[i,:]) for i in barlabels], index=barlabels)
    norm_gdfs = age_normdf
    cum_gdfs = age_cumdf
    for i, col in enumerate(age_df.columns):
        height = norm_gdfs[col]
        starts = cum_gdfs[col] - height
        rects = ax.barh(barlabels, height, left=starts, height=0.9, color=colors[i], edgecolor='white', linewidth=0.5,
                        label=col, alpha=1)
    labels = [f'{val:.2f}' for val in age_norm1df.loc[:, age_norm1df.columns[0]]]
    if bar_label == 'first':
      counts = [int(val) for val in age_normdf.loc[:, age_norm1df.columns[0]]]
    else:
      counts = [int(val) for val in all_df]
    ax.bar_label(rects, counts, 
                 label_type='edge', color='#7f3730', fontsize=10)
        # ax.bar_label(rects,age_df.loc[:, col], label_type='center', color='lightgrey', fontsize=14)
    ax.legend( bbox_to_anchor=(1, 0.3), 
              handletextpad=0.5, frameon=False,
                          borderpad=0.4,
                          columnspacing=1,
                          handlelength=0.65,
              loc='lower left')
    ax.set_ylabel(xlabel)
    def ks(x, pos):
        return '%1.1fk' % (x*1e-3)
    ax.xaxis.set_major_formatter(FuncFormatter(ks))
    # ax.set_yticks(barlabels)
    # ax.set_xticklabels(labels=barlabels,rotation=xrotate)
    ax.spines.top.set_visible(False)
    ax.spines.right.set_visible(False)

def build_EBV_gtf():
  """
  生成EBV的gtf文件。

  """
  EBV_gtf = pd.read_csv("~/yard/refdata-cellranger-EBV/reference_sources/chrEBV_Akata_inverted_refined_genes_plus_features_annotation_cleaned_2.gtf", sep='\t',header=None, comment='#', names=['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
  GRC_gtf = pd.read_csv("~/yard/refdata-cellranger-GRCh38-EBV-2024-A/reference_resources/genes.gtf", sep='\t',header=None, comment='#', names=['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
  merged_gtf = pd.concat([GRC_gtf, EBV_gtf], ignore_index=True)
  merged_gtf.to_csv('~/yard/refdata-cellranger-GRCh38-EBV-2024-A/reference_resources/merged.gtf', sep='\t', header=False, index=False)
  recordsA = list(SeqIO.parse("/home/rzh/yard/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa", "fasta"))
  recordsB = list(SeqIO.parse("/home/rzh/yard/refdata-cellranger-EBV/reference_sources/chrEBV_Akata_inverted_2.fa", "fasta"))
  merged_records = recordsA + recordsB
  with open('/home/rzh/yard/refdata-cellranger-GRCh38-EBV-2024-A/reference_resources/merged_file.fa', 'w') as outfile:
    SeqIO.write(merged_records, outfile, "fasta")
