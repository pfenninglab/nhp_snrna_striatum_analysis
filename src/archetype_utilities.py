import scanpy as sc
import pandas as pd
import numpy as np
from scipy.io import mmread
import anndata
from VLAD.vlad import vlad, vlad_centers_given_2
from rdd import rdd
from skimage.filters import threshold_otsu
from sklearn.manifold import SpectralEmbedding
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr, pearsonr
import seaborn as sns
from statsmodels.stats.multitest import multipletests
import warnings
from sklearn.cluster import AgglomerativeClustering

def get_count_matrix(monkey_a_annot,monkey_a_mat,cell_types_chosen,gene_indices):
    monkey_a_indices = np.logical_or(monkey_a_annot["cell_type"].values == cell_types_chosen[0],
                                        monkey_a_annot["cell_type"].values == cell_types_chosen[1])
    monkey_a_dense = np.asarray(monkey_a_mat.todense())[monkey_a_indices,:][:,gene_indices]
    return monkey_a_dense


#hot start archetype projection based on already selected cell types
def archetype_projection(monkey_a_mat, monkey_b_mat, monkey_a_annot, monkey_b_annot,
                         cell_types_chosen, gene_names_vec):
    monkey_a_indices = np.logical_or(monkey_a_annot["cell_type"].values == cell_types_chosen[0],
                                     monkey_a_annot["cell_type"].values == cell_types_chosen[1])
    monkey_a_dense = np.asarray(monkey_a_mat.todense())[monkey_a_indices, :]
    monkey_a_ann = anndata.AnnData(monkey_a_dense, obs=monkey_a_annot.loc[monkey_a_indices, :])
    monkey_a_ann.var_names = gene_names_vec
    sc.pp.normalize_total(monkey_a_ann)
    sc.pp.pca(monkey_a_ann)
    sc.pp.neighbors(monkey_a_ann)
    sc.tl.umap(monkey_a_ann)

    a = (monkey_a_ann.obsp['connectivities'] + monkey_a_ann.obsp['connectivities'].transpose()) / 2
    spec = SpectralEmbedding(affinity="precomputed")
    spectral_a = spec.fit_transform(a)
    kmeans = KMeans(n_clusters=2)
    kmeans.fit(spectral_a)
    clusassign = kmeans.fit_predict(spectral_a)
    min_dist = np.min(cdist(spectral_a, kmeans.cluster_centers_, 'euclidean'), axis=1)
    Y = pd.DataFrame(min_dist, columns=['Center_euclidean_dist'])
    Z = pd.DataFrame(clusassign, columns=['cluster_ID'])
    PAP = pd.concat([Y, Z], axis=1)
    grouped = PAP.groupby(['cluster_ID'])
    grouped.idxmin()
    centers = np.zeros(a.shape[0])
    centers[grouped.idxmin().values[0][0]] = 1
    centers[grouped.idxmin().values[1][0]] = 1
    monkey_a_ann.obs["centers"] = pd.Categorical(centers)
    center_points = (grouped.idxmin().values[0][0], grouped.idxmin().values[1][0])
    cell_labels = np.asarray(monkey_a_ann.obs["cell_type"].values)
    # cell_labels[monkey_a_ann.obs["centers"].values.astype(bool)] = "center"
    # monkey_a_ann.obs["label"] = pd.Categorical(cell_labels)
    # sc.pl.umap(monkey_a_ann,color="label")

    monkey_a_beta = vlad_centers_given_2(X=monkey_a_dense, K=2, center_indices=center_points, alpha=.97, pois=True)
    monkey_a_beta_inv = np.linalg.pinv(monkey_a_beta)
    monkey_a_weights = monkey_a_ann.X.dot(monkey_a_beta_inv)
    monkey_a_ann.obs["weight_1"] = monkey_a_weights[:, 0]
    monkey_a_ann.obs["weight_2"] = monkey_a_weights[:, 1]
    # sc.pl.scatter(monkey_a_ann,"weight_1","weight_2",color="cell_type")

    monkey_b_indices = np.logical_or(monkey_b_annot["cell_type"].values == cell_types_chosen[0],
                                     monkey_b_annot["cell_type"].values == cell_types_chosen[1])
    monkey_b_dense = np.asarray(monkey_b_mat.todense())[monkey_b_indices, :]
    monkey_b_ann = anndata.AnnData(monkey_b_dense, obs=monkey_b_annot.loc[monkey_b_indices, :])
    monkey_b_ann.var_names = gene_names_vec
    # sc.pp.normalize_total(monkey_b_ann)
    # sc.pp.pca(monkey_b_ann)
    # sc.pp.neighbors(monkey_b_ann)
    # sc.tl.umap(monkey_b_ann)
    monkey_b_weights = monkey_b_ann.X.dot(monkey_a_beta_inv)
    monkey_b_ann.obs["weight_1"] = monkey_b_weights[:, 0]
    monkey_b_ann.obs["weight_2"] = monkey_b_weights[:, 1]
    sc.pl.scatter(monkey_b_ann, "weight_1", "weight_2", color="cell_type")

    return monkey_a_ann, monkey_b_ann, monkey_a_beta, monkey_a_beta_inv

#unsupervised projection withk archetypes
def archetype_multi_projection(monkey_a_mat, monkey_b_mat, monkey_a_annot, monkey_b_annot,
                               cell_types_chosen, gene_names_vec, num_archetypes):
    monkey_a_indices = np.logical_or(monkey_a_annot["cell_type"].values == cell_types_chosen[0],
                                     monkey_a_annot["cell_type"].values == cell_types_chosen[1])
    monkey_a_dense = np.asarray(monkey_a_mat.todense())[monkey_a_indices, :]
    monkey_a_ann = anndata.AnnData(monkey_a_dense, obs=monkey_a_annot.loc[monkey_a_indices, :])
    monkey_a_ann.var_names = gene_names_vec
    sc.pp.normalize_total(monkey_a_ann)
    sc.pp.pca(monkey_a_ann)
    sc.pp.neighbors(monkey_a_ann)
    sc.tl.umap(monkey_a_ann)
    sc.pl.umap(monkey_a_ann, color="cell_type")

    monkey_a_beta = vlad(X=monkey_a_dense, K=num_archetypes, alpha=.97, pois=True)
    monkey_a_beta_inv = np.linalg.pinv(monkey_a_beta)
    monkey_a_weights = monkey_a_ann.X.dot(monkey_a_beta_inv)
    for i in range(num_archetypes):
        monkey_a_ann.obs["weight_" + str(i + 1)] = monkey_a_weights[:, i]

    monkey_b_indices = np.logical_or(monkey_b_annot["cell_type"].values == cell_types_chosen[0],
                                     monkey_b_annot["cell_type"].values == cell_types_chosen[1])
    monkey_b_dense = np.asarray(monkey_b_mat.todense())[monkey_b_indices, :]
    monkey_b_ann = anndata.AnnData(monkey_b_dense, obs=monkey_b_annot.loc[monkey_b_indices, :])
    monkey_b_ann.var_names = gene_names_vec
    sc.pp.normalize_total(monkey_b_ann)
    sc.pp.pca(monkey_b_ann)
    sc.pp.neighbors(monkey_b_ann)
    sc.tl.umap(monkey_b_ann)
    monkey_b_weights = monkey_b_ann.X.dot(monkey_a_beta_inv)
    for i in range(num_archetypes):
        monkey_b_ann.obs["weight_" + str(i + 1)] = monkey_b_weights[:, i]

    return monkey_a_ann, monkey_b_ann, monkey_a_beta, monkey_a_beta_inv


def vlad_mstd(data_mat,k_range,num_repeats,alpha=.97):
    stability_list = []
    for k in k_range:
        print(k)
        comp_list = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for i in range(num_repeats):
                beta = vlad(X=data_mat, K=k,alpha=alpha,pois=True)
                comp_list.append(beta)
        comp_mat = np.concatenate(comp_list,axis=0)
        comp_abs_corr = np.nan_to_num(np.corrcoef(comp_mat), posinf=1, neginf=-1)
        #comp_abs_corr = np.abs(cosine_similarity(comp_mat))
        agg_clust = AgglomerativeClustering(affinity="precomputed",n_clusters=k,linkage="average")
        #agg_clust = SpectralClustering(affinity="precomputed",n_clusters=k)
        agg_clust.fit(1.0-comp_abs_corr)
        print(agg_clust.labels_)
        s_index =  stability_index(agg_clust.labels_,comp_abs_corr,k)
        stability_list.append(s_index)
    mstd_df = pd.DataFrame({"n_comps":list(k_range),"stability":stability_list})
    return mstd_df



def vlad_rep(data_mat,k_range,num_repeats,alpha):
    stability_dict = {}
    for k in k_range:
        print(k)
        comp_list = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for i in range(num_repeats):
                beta = vlad(X=data_mat, K=k,alpha=alpha,pois=True)
                comp_list.append(beta)
        stability_dict[k] = comp_list
    return stability_dict



def stability_index(clustering, abs_dist_mat, num_clusts):
    st_index_sum = 0
    num_comps = abs_dist_mat.shape[0]
    clust_values = np.unique(clustering)
    for k in clust_values:
        pos_part = 0
        neg_part = 0
        for i in range(num_comps):
            for j in range(i + j, num_comps):
                if clustering[i] == clustering[j]:
                    pos_part += abs_dist_mat[i, j]
                else:
                    neg_part += abs_dist_mat[i, j]
        clust_size = np.sum(clustering == k)
        # print(clust_size)
        pos_coeff = 1.0 / np.square(clust_size)
        neg_coeff = 1.0 / (clust_size * (num_comps - clust_size))
        st_index_sum += (pos_coeff * pos_part) - (neg_coeff * neg_part)

    stability_index = (1.0 / num_clusts) * st_index_sum
    return stability_index
