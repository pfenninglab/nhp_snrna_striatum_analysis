import matplotlib.pyplot as plt
import seaborn
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.io import mmread
import anndata
from VLAD.vlad import vlad
import pickle
from skimage.filters import threshold_otsu
from scipy.stats import pointbiserialr, spearmanr, pearsonr
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import statsmodels.formula.api as smf
import statsmodels.api as sm

def dense_hist(x, **kwargs):
    sns.histplot(x, stat="probability", **kwargs).set_title('pair')


def calc_rdd_values_glm_permute(data_df, gene_list, weight_vec, posmark, ngenes_vec):
    pvals = []
    tvals = []
    remove_list = []
    gene_list_fixed = []
    for gene in gene_list:
        np.random.shuffle(posmark)
        test_df = pd.DataFrame(
            {"weight_vec": weight_vec, "posmark": posmark, "ngenes": ngenes_vec, "gene": data_df[gene].values})
        try:
            rdd_fit = smf.glm(formula='gene ~ posmark+weight_vec + ngenes', data=test_df,
                              family=sm.families.Poisson()).fit()
            # _,pvals_corrected,_,_ = multipletests(rdd_fit.pvalues[1] ,method="fdr_bh")
            pvals.append(rdd_fit.pvalues[1])
            tvals.append(rdd_fit.tvalues[1])
            gene_list_fixed.append(gene)
        except:
            pvals.append(1.0)
            tvals.append(0.0)
            gene_list_fixed.append(gene)
            print("regression issue")

    _, pvals_corrected, _, _ = multipletests(pvals, method="fdr_bh")
    rdd_df = pd.DataFrame(
        {"gene": gene_list_fixed, "p_value": pvals, "p_value_corrected": pvals_corrected, "t_value": tvals})
    return rdd_df


def corr_rdd_df_glm_permute(monkey_a_ann, monkey_b_ann, monkey_b_annot, monkey_b_mat, archetype_label, gene_names_vec,
                            iteration, pair_key):
    cell_types_chosen = np.unique(monkey_a_ann.obs["cell_type"].values)

    corr_pval_calc = lambda x: pearsonr(x, monkey_a_ann.obs[archetype_label].values)
    archetype_corr = np.apply_along_axis(corr_pval_calc, 0, monkey_a_ann.X)
    archetype_corr_series = pd.Series(archetype_corr[1, :], index=gene_names_vec)
    archetype_corr_series = archetype_corr_series.fillna(1)
    _, archetype_corr_pvals_corrected, _, _ = multipletests(archetype_corr_series.values, method="fdr_bh")
    archetype_corr_corrected_series = pd.Series(archetype_corr_pvals_corrected, index=gene_names_vec)
    genes_selected = archetype_corr_corrected_series[archetype_corr_pvals_corrected < 0.05].index.values
    # print(genes_selected)
    # archetype_corr_values = archetype_corr[0,:][archetype_corr_pvals_corrected<0.05]
    corr_values_series = pd.Series(archetype_corr[0, :], index=gene_names_vec)

    monkey_b_smat = get_count_matrix(monkey_b_annot, monkey_b_mat, cell_types_chosen,
                                     archetype_corr_pvals_corrected < 0.05)
    ngenes_vec = np.sum(monkey_b_smat > 0, axis=1)
    monkey_b_df = pd.DataFrame(monkey_b_smat, columns=genes_selected)

    # print("!!!!Threshold:!!!!!")
    thresh = threshold_otsu(monkey_b_ann.obs[archetype_label].values)
    wvec = monkey_b_ann.obs[archetype_label].values - thresh
    posthresh = wvec > 0
    # np.random.shuffle(posthresh)
    # print(thresh)
    rdd_df = calc_rdd_values_glm_permute(monkey_b_df, genes_selected, wvec, posthresh, ngenes_vec)
    rdd_df["pearson_pvalue"] = archetype_corr_corrected_series[rdd_df["gene"].values].values
    rdd_df["pearson_corr"] = corr_values_series[rdd_df["gene"].values].values
    rdd_df["neg_log_pvalue"] = -1 * np.log(rdd_df["p_value"])
    rdd_df["neg_log_pvalue_corrected"] = -1 * np.log(rdd_df["p_value_corrected"])
    rdd_df["neg_log_spearman_pvalue"] = -1 * np.log(rdd_df["pearson_pvalue"])
    rdd_df["iteration"] = iteration
    rdd_df["pair"] = pair_key
    sns.histplot(rdd_df["p_value"].values)
    sns.scatterplot(data=rdd_df, x="neg_log_spearman_pvalue", y="neg_log_pvalue")
    return rdd_df
