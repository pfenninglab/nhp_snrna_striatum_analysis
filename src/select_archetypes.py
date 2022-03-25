import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial.distance import hamming
from scipy.optimize import linear_sum_assignment

from scipy.cluster.hierarchy import linkage, cophenet
from scipy.spatial.distance import pdist, squareform, cdist
from scipy.optimize import linear_sum_assignment

#adopted from https://www.pnas.org/doi/10.1073/pnas.0308531101
def cophenetic(A):
    avec = np.array([A[i, j] for i in range(A.shape[0] - 1)
                     for j in range(i + 1, A.shape[1])])
    # consensus entries are similarities, conversion to distances
    Y = 1 - avec
    Z = linkage(Y, method='average')
    return cophenet(Z, Y)[0]

#adopted from https://academic.oup.com/bioinformatics/article/23/12/1495/225472
def dispersion(A):
    return np.sum(4 * np.square(A - .5)) / A.size

#adopted from https://www.pnas.org/doi/10.1073/pnas.0308531101
def consensus_matrix_binary(beta_list):
    max_list = []
    for beta_i in beta_list:
        max_list.append(np.argmax(beta_i, axis=0))
    max_mat = np.asarray(max_list)
    return 1.0 - squareform(pdist(max_mat, hamming))

#adopted from https://www.pnas.org/doi/10.1073/pnas.0308531101
def consensus_matrix_continous(beta_list, metric="cosine"):
    avg_mat = None
    for beta_i in beta_list:
        mat_i = 1.0 - squareform(pdist(beta_list, metric=metric))
        if avg_mat is None:
            avg_mat = mat_i
        else:
            avg_mat = avg_mat + mat_i
    return avg_mat / len(beta_list)

#Uses the hungararian algorithm, this won't scale for high N but is fine with low N
def maximum_matching(beta_list, metric="cosine"):
    run_sum = 0
    num_betas = len(beta_list)
    for i in range(num_betas):
        for j in range(i + 1, num_betas):
            dist_mat = cdist(beta_list[i], beta_list[j], metric=metric)
            min_dist_idx = linear_sum_assignment(dist_mat)
            score = np.sum(np.square(dist_mat[min_dist_idx[0], min_dist_idx[1]])) / dist_mat.shape[0]
            run_sum += score
    return run_sum / (np.square(num_betas))
