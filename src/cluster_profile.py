from scipy.stats import hypergeom
import pandas as pd
import numpy as np
from anndata import AnnData
import scanpy as sc
from typing import Dict
from IPython.display import display, HTML




def profile_anndata(ann:AnnData,cell_type_marker_dict:Dict,cluster_key="leiden",pval_cutoff=0.1,
                    print_to_notebook=True):
    '''
    Function profiles single nuclei rna-seq clustering based on a marker gene dictionary
    A hypergeometric test use to perform enrichments
    :param ann: AnnData containing single cell gene expression
    :param cell_type_marker_dict: dict with keys as cell type names and values as sets of genes
    :param cluster_key: key of clustering index in ann
    :param pval_cutoff: only consider enrichments with pvalues smaller than this cutoff
    :param print_to_notebook: whether the results are being print to a jupyter notebook
    :return: dict with keys as clusters and values as the assigned cell type by most significant enrichment
    '''
    cluster_assign_dict = {}
    all_gene_set =set(ann.var_names.tolist())
    #bonferroni correction is based on the total number of tests
    bonferroni_correction = np.unique(ann.obs[cluster_key].values).size * len(cell_type_marker_dict.keys())
    for cluster in np.unique(ann.obs[cluster_key].values):
        print("------------------------------------------------------------------------------")

        to_print = "CLUSTER: " + str(cluster)
        print(to_print)
        # calculate marker genes for each cluster
        is_mark = [ctype == cluster for ctype in ann.obs[cluster_key].values]
        ann.obs["mark"] = pd.Categorical(
            values=is_mark,
            categories=[True, False])
        sc.tl.rank_genes_groups(ann, "mark")
        gene_set = set([tup[0] for tup in ann.uns["rank_genes_groups"]["names"]])

        best_cell_type = "None"
        best_pvalue = pval_cutoff
        pvals = []
        cell_types = []
        genes_found = []
        #for each cell type in the given marker gene dictionary perform an intersection of marker genes with the
        #cluster and then calculate statistical significance with the hpyer geometric test
        for cell_type in cell_type_marker_dict:
            # intersect with cluster markers
            inter = gene_set.intersection(cell_type_marker_dict[cell_type])
            #intersection with background
            all_inter = all_gene_set.intersection(cell_type_marker_dict[cell_type])
            pval = hypergeom.sf(len(inter) - 1, len(all_gene_set), len(all_inter),
                                len(gene_set)) * bonferroni_correction
            if pval < pval_cutoff and len(inter) > 0:
                pvals.append(pval)
                cell_types.append(cell_type)
                genes_found.append(','.join(inter))
            if pval < best_pvalue:
                best_pvalue = pval,
                best_cell_type = cell_type
        #display enrichments
        result_dict = {"p_value": pvals, "label": cell_types, "query_match": genes_found}
        if print_to_notebook:  display(HTML(pd.DataFrame(result_dict).to_html()))
        else: print(pd.DataFrame(result_dict))
        result_assignment = "Cluster " + str(cluster) + " has been labeled as " + best_cell_type
        print(result_assignment)
        cluster_assign_dict[cluster] = best_cell_type
    return cluster_assign_dict



