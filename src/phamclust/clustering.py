from sklearn.cluster import AgglomerativeClustering


def hierarchical(matrix, linkage, eps=None, n_clusters=None):
    """Perform hierarchical clustering of the matrix nodes.

    Returns cluster sub-matrices in descending size order.

    NOTE: parameters `eps` and `n_clusters` are mutually exclusive.

    :param matrix: the matrix whose nodes are to be clustered
    :type matrix: phamclust.matrix.SymMatrix
    :param linkage: the kind of linkage to use for clustering
    :type linkage: str
    :param eps: the distance threshold to use for clustering
    :type eps: float
    :param n_clusters: the number of clusters to return
    :type n_clusters: int
    :return: clusters
    """
    if len(matrix) == 1:
        return [matrix]

    # Export distance matrix
    if not matrix.is_distance:
        raise ValueError(f"matrix must be a distance matrix")
    dist_mat = matrix.to_ndarray()

    if eps is None and not n_clusters:
        raise ValueError(f"need either threshold or n_clusters to proceed")
    elif eps and n_clusters:
        raise ValueError(f"threshold and n_clusters are mutually exclusive")

    # Build the clustering object using indicated eps/n_clusters
    c = AgglomerativeClustering(metric="precomputed", linkage=linkage,
                                distance_threshold=eps, n_clusters=n_clusters)

    temp_clusters = dict()
    for node_name, label in zip(matrix.nodes, c.fit_predict(dist_mat)):
        if label in temp_clusters:
            temp_clusters[label].add(node_name)
        else:
            temp_clusters[label] = {node_name}

    # sklearn API gives singletons labels just like any other cluster
    clusters = list()
    for node_names in temp_clusters.values():
        clusters.append(matrix.extract_submatrix(list(node_names)))

    # Return clusters in reverse order of size (largest to smallest)
    return sorted(clusters, reverse=True)
