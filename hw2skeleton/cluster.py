from hw2skeleton.utils import Atom, Residue, ActiveSite


def extract_aa_seq(active_site):
    return [letters[x.type[0:3]] for x in active_site.residues]


def compute_similarity(site_a, site_b):
    """
    Caculate the similarity base solely on AA sequence similarity
    Divide the number of common residues (with multiplicity)
    by the greater number of each type of residue present.
    i.e., s('ABC', 'ADA') =
    A - 1 in common, 2 max
    B - 0 in common, 1 max
    C - 0 in common, 1 max
    D - 0 in common, 1 max
    (1 + 0 + 0 + 0) / (2 + 1 + 1 + 1) = 1 / 5 = 0.2
    """
    seq1 = extract_aa_seq(site_a)
    seq2 = extract_aa_seq(site_b)
    seen_res = set(seq1 + seq2)
    total = 0
    common = 0
    for res in seen_res:
        in_1 = seq1.count(res)
        in_2 = seq2.count(res)
        if in_1 < in_2:
            common += in_1
            total += in_2
        else:
            common += in_2
            total += in_1
    similarity = common / total
    return similarity


def cluster_by_partitioning(active_sites):
"""
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    Algorithm: This is an implementation of the affinity propagation algorithm,
               as described in 'Clustering by Passing Messages between Data Points',
               Frey and Dueck, Science  16 Feb 2007: Vol. 315, Issue 5814, pp. 972-976
               DOI: 10.1126/science.1136800
               
               Clusters are based on so-called 'exemplary points', which maximize the
               availability and responsibility of points in their clusters. 
               
               Here, responsibility can be understood as the likelihood that one point 
               should choose a second point to be exemplary (cluster center) over all other
               candidate exemplars. 
               (likelihood exemplar(a) = b given all other possible exemplars)
               
               Availability is the likelihood that one point should be chosen by a second
               point to be its exemplar, given all of its existing connections.
               (likelihood exemplar(b) = a given all other members of cluster a)
    """
    sim_mat = np.asarray([np.asarray([compute_similarity(x, y) for x in active_sites]) for y in active_sites])
    a = np.zeros((sim_mat.shape[0], sim_mat.shape[1]))
    r = np.zeros((sim_mat.shape[0], sim_mat.shape[1]))
    damp = 0.4
    sim_mat += np.random.rand(sim_mat.shape[0], sim_mat.shape[1]
                              ) * 1.e-12 * (np.amax(sim_mat) - np.amin(sim_mat))
    for i in range(0, 300):
        r_old = r
        a_s = a + sim_mat
        y, ind = np.max(a_s, axis=1), np.argmax(a_s, axis=1)
        for j in range(sim_mat.shape[0]):
            a_s[j, ind[j]] = -1 * np.inf
        y_2, ind_2 = np.max(a_s, axis=1), np.argmax(a_s, axis=1)
        r = sim_mat - np.tile(y, (sim_mat.shape[0], 1))
        for j in range(sim_mat.shape[0]):
            r[j][ind[j]] = sim_mat[j][ind[j]] - y_2[j]
        r = (1 - damp) * r + damp * r_old
        
        a_old = a
        r_p = r
        r_p[r_p < 0] = 0
        for j in range(sim_mat.shape[0]):
            r_p[j][j] = r[j][j]
        a = np.tile(np.sum(r_p, axis=0), (sim_mat.shape[0], 1)) - r_p
        da = np.diagonal(a)
        a[a > 0] = 0
        for j in range(sim_mat.shape[0]):
            a[j][j] = da[j]
        a = (1 - damp) * a + damp * a_old
        
    e = r + a
    ind = np.where(np.diagonal(e) != 0)[0]
    k = ind.size
    c = np.argmax(sim_mat[:, ind], axis=1)
    c[ind] = np.arange(k)
    for j in range(k):
        tmp = np.where(c == j)[0]
        dest = np.argmax(np.sum(sim_mat[tmp[:, np.newaxis], tmp], axis=0))
        ind[j] = tmp[dest]
    c = np.argmax(sim_mat[:, ind], axis=1)
    c[ind] = np.arange(k)
    labels = ind[c]
    clust_cent_ind = np.unique(labels)
    labels = np.searchsorted(clust_cent_ind, labels)
    as_clusters = [(labels[i], active_sites[i].name) for i in range(
        0, len(active_sites))]
    as_clusters = sorted(as_clusters)
    clusters = [[]]
    current_cluster = 0
    names = [x.name for x in active_sites]
    for i in range(0, len(active_sites)):
        if as_clusters[i][0] != current_cluster:
            clusters.append([])
        clusters[-1].append(active_sites[names.index(as_clusters[i][1])])
        current_cluster = as_clusters[i][0]
        return clusters


# this function flattens a nested list
# recursion is not used to avoid the recursion barrier
def flatten(l):
    if isinstance(l, list):
        go_deeper = True
        while go_deeper:
            tmp = []
            go_deeper = False
            for x in l:
                if isinstance(x, list):
                    tmp += x
                    go_deeper = True
                else:
                    tmp.append(x)
            l = tmp
    return l


def cluster_hierarchically(active_sites):
"""
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """
    dist_mat = 1.0 - np.asarray([np.asarray([compute_similarity(x, y) for x in active_sites]
                                    ) for y in active_sites])
    in_which_cluster = np.array(list(range(0, dist_mat.shape[0])))
    clusters = list(range(0, dist_mat.shape[0]))
    for i in range(0, dist_mat.shape[0]):
        dist_mat[i][i] = np.inf
    while(len(clusters)) > 2:
        nearest_neighbor = np.argmin(dist_mat)
        nn_i, nn_j = int(nearest_neighbor / dist_mat.shape[0]), nearest_neighbor % dist_mat.shape[0]
        if in_which_cluster[nn_i] != in_which_cluster[nn_j]:
            # join the clusters
            old_clust = clusters[in_which_cluster[nn_j]]
            old_nnj = in_which_cluster[nn_j]
            clusters[in_which_cluster[nn_i]] = [
                clusters[in_which_cluster[nn_j]], clusters[in_which_cluster[nn_i]]]
            # set the location of the new component's members as the old component's location
            in_which_cluster[flatten(clusters[in_which_cluster[nn_j]])] = in_which_cluster[nn_i]
            # remove the old, independent cluster
            clusters.remove(old_clust)
            # when a new cluster is created, shift the index of all clusters 
            # to the right of the greater position by -1
            for i in range(0, len(in_which_cluster)):
                if in_which_cluster[i] > old_nnj:
                    in_which_cluster[i] -= 1
        dist_mat[nn_i][nn_j] = np.inf
        dist_mat[nn_j][nn_i] = np.inf
    return clusters


# quality is defined as the average over all combinations of clusters
# of the product of the geometric mean of the intra-cluster similarities
# and the inter-cluster dissimilarity (1 - similarity)
def clustering_quality(list_of_clusters, active_sites):
    num_to_ind = {active_sites[i]: i for i in range(0, len(active_sites))}
    sim_mat = np.asarray([np.asarray([compute_similarity(x, y) for x in active_sites]
                                ) for y in active_sites])
    n = len(list_of_clusters)
    intra_sim = [(np.sum([[sim_mat[num_to_ind[u]][num_to_ind[v]] for u in x] for v in x
                         ])) / (len(x) ** 2) for x in list_of_clusters]
    quality = 0
    for i in range(0, len(list_of_clusters)):
        for j in range(i, len(list_of_clusters)):
            quality += np.sqrt(intra_sim[i] * intra_sim[j]) * (
                1.0 - np.sum([[sim_mat[num_to_ind[u]][num_to_ind[v]
                            ] for u in list_of_clusters[i]] for v in list_of_clusters[j]]
                          ) / (len(list_of_clusters[i]) + len(list_of_clusters[j])))
    quality = (1 / (n * (n - 1))) * quality
    return quality
