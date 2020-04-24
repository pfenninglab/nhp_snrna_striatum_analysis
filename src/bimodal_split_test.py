import scanpy as sc
import pandas as pd
import numpy as np
from tableone import modality
from statsmodels.stats.multitest import multipletests
import argparse
import anndata

def diptest_dip(vec:np.ndarray):
    '''
    Wrapper over hartigan's dip test
    :param vec: vector of real values
    :return: Hartigan's dip test p-value
    '''
    return modality.hartigan_diptest(vec)



def diptest_dip_zinf(vec:np.ndarray):
    '''
    Outputs hartigan's diptest assuming a zero inflated vector
    :param vec: vector of real values
    :return: Hartigan's dip test p-value on non-zero values
    '''
    if np.sum(vec!=0)<10:
        return 1.0
    return modality.hartigan_diptest(vec[vec!=0])

def choose_pvals(pval_vec:np.ndarray,zinf_pval_vec:np.ndarray):
    '''
    method to decide whether to use the normal or zero-inflated statistics based on pvalues
    :param pval_vec: vector of gene dip test pvalues
    :param zinf_pval_vec:  vector of gene dip test pvalues assuming zero inflation
    :return: vector of uncorrected gene pvalues
    '''
    pvals = np.zeros(len(pval_vec))
    for i in range(len(pval_vec)):
        if zinf_pval_vec[i]>pval_vec[i]*2:
            pvals[i] = -1*zinf_pval_vec[i]
        else: pvals[i] = pval_vec[i]
    return pvals

def combine_tests(vec:np.ndarray):
    '''
    Helper method to enforce the following logic, if a gene is zero-inflated in both splits
    then the full vector should use the zero inflated pvalue too, otherwise it should use the regular pvalue
    the full vector split pvalue should only be reported if both splits are unimodal

    '''
    if np.abs(vec[0])>.1 and np.abs(vec[1])>.1:
        if vec[0]<0 and vec[1]<0:
            return np.abs(vec[3])
        else:
            return vec[2]
    else: return 1.0



def split_pvals(mat:np.ndarray,split_vec:np.ndarray):
    '''

    :param mat: single nuclei gene expression matrix, formated as : nuclei x genes
    :param split_vec: boolean vector representing the cluster split
    :return: Benjamini Hochberg corrected pvalues for each gene for the bimodal split test
    '''
    mat_a = mat[split_vec,:]
    mat_b = mat[np.logical_not(split_vec)]

    dip_a = np.apply_along_axis(diptest_dip, 0, mat_a)
    zinf_a = np.apply_along_axis(diptest_dip_zinf, 0, mat_a)

    dip_a_fixed = choose_pvals(dip_a,zinf_a)


    dip_b = np.apply_along_axis(diptest_dip, 0, mat_b)

    zinf_b = np.apply_along_axis(diptest_dip_zinf, 0, mat_b)
    dip_b_fixed = choose_pvals(dip_b, zinf_b)


    dip_split = np.apply_along_axis(diptest_dip,0,mat)
    #_, corrected_dip_split, _, _ = multipletests(dip_split, method="fdr_bh")
    dip_split_zinf = np.apply_along_axis(diptest_dip_zinf,0,mat)

    #_, corrected_dip_split_zinf, _, _ = multipletests(dip_split_zinf, method="fdr_bh")
    pval_arr = np.vstack([dip_a_fixed,dip_b_fixed,dip_split,dip_split_zinf])
    raw_pvals = np.apply_along_axis(combine_tests,0,pval_arr)

    _,corrected_pval_arr,_,_ = multipletests(raw_pvals,method="fdr_bh")

    return corrected_pval_arr

def split_pvals_ann(ann:anndata.AnnData,split_vec):
    '''

    :param ann: Annotated dataframe of single nuclei gene expression values
    :param split_vec: boolean vector representing the cluster split
    :return: pandas Series with genes as the index and values of Benjamini Hochberg corrected
     pvalues for each gene for the bimodal split test
    '''
    sc.pp.filter_genes(ann, min_cells=10)
    pvals = pd.DataFrame(ann.X, columns=ann.var_names)
    return pd.Series(pvals,index=ann.var_names)


def bootstrap_analysis(ann:anndata.AnnData,num_bootstraps,outfile):
    '''
    Performs bootstrap analysis to generate a null distribution
    :param ann: Annotated dataframe of single nuclei gene expression values
    :param num_bootstraps: Number of bootstrap to perform, the higher the number the
    :param outfile: outfile to write values for null distribution
    :return: None
    '''
    sc.pp.filter_genes(ann, min_cells=10)
    exp_df = pd.DataFrame(ann.X, columns=ann.var_names)

    for i in range(num_bootstraps):
        if i%100==0:
            print(i)
        random_select = np.random.rand(exp_df.shape[0])>.5
        random_pvals = split_pvals(exp_df.values, random_select)
        line = str(np.sum(random_pvals<0.05)) +"\n"
        if i<10:
            print(line)
        f = open(outfile, "a+")
        f.write(line)
        f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="Input h5ad of AnnData file path")
    parser.add_argument("-m", help="Mode: should be either b for boostrap or s for split results",default='s')
    parser.add_argument("-n", help="Number of bootstraps ", type=int, default=10000)
    parser.add_argument("-v", help="Object in AnnData that represents the cluster split ", type=int, default=10000)
    parser.add_argument("-o", help="Outfile path")

    args = parser.parse_args()

    mode = pd.read_csv(args.m)
    ann =  anndata.read_h5ad(args.i)

    if mode=='s':
        split_vec =  ann.obs[args.v].values
        gene_pvalues = split_pvals_ann(ann,split_vec)
        pd.DataFrame(gene_pvalues).to_csv(args.o)
    elif mode=='b':
        num_bootstraps = args.n
        bootstrap_analysis(ann,num_bootstraps,args.o)
    else: print("Mode not supported: should be either b for boostrap analysis or s for split results")









