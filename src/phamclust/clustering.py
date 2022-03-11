"""Functions for clustering the genomes in a symmetric matrix."""

import queue

from phamclust.SymMatrix import SymMatrix


def single_linkage(matrix, eps=None):
    """Cluster matrix with single-linkage to the indicated threshold.

    Implementation uses dbscan with minpts=1.

    :param matrix: the matrix whose nodes need clustering
    :type matrix: SymMatrix
    :param eps: similarity threshold for clustering
    :type eps: float
    :return: clusters
    """
    return dbscan(matrix, eps, minpts=1)


def hierarchical(matrix, eps=None):
    """Cluster matrix hierarchically to the indicated threshold.

    Implementation uses dbscan with minpts=2.

    :param matrix: the matrix whose nodes need clustering
    :type matrix: SymMatrix
    :param eps: similarity threshold for clustering
    :type eps: float
    :return: clusters
    """
    return dbscan(matrix, eps, minpts=2)


def dbscan(matrix, eps=None, minpts=3):
    """Cluster matrix with dbscan to the indicated threshold.

    :param matrix: the matrix whose nodes need clustering
    :type matrix: SymMatrix
    :param eps: similarity threshold for clustering
    :type eps: float
    :param minpts: minimum number of neighbors to build a cluster
    :type minpts: int
    :return: clusters
    """
    if eps is None:
        eps = matrix.epsilon

    counter = 0
    cluster_lookup = dict()
    for source_node in matrix.node_names:
        # Skip nodes we clustered already
        if cluster_lookup.get(source_node, -1) != -1:
            continue

        neighbors = matrix.nearest_neighbors(source_node, eps)
        if len(neighbors) < minpts:
            cluster_lookup[source_node] = None
            continue

        counter += 1
        cluster_lookup[source_node] = counter

        seeds = queue.Queue()
        for neighbor in neighbors:
            seeds.put(neighbor)

        while not seeds.empty():
            target = seeds.get()

            target_cluster = cluster_lookup.get(target, -1)
            if not target_cluster:
                cluster_lookup[target] = counter
            if target_cluster != -1:
                continue
            cluster_lookup[target] = counter
            neighbors = matrix.nearest_neighbors(target, eps)

            if len(neighbors) >= minpts:
                for neighbor in neighbors:
                    seeds.put(neighbor)

    clusters = dict()
    for node_name, cluster in cluster_lookup.items():
        if not cluster:
            counter += 1
            clusters[counter] = [node_name]
        else:
            cluster_nodes = clusters.get(cluster, list())
            cluster_nodes.append(node_name)
            clusters[cluster] = cluster_nodes
    for cluster, cluster_nodes in clusters.items():
        cluster_matrix = matrix.extract_submatrix(cluster_nodes)
        edges = cluster_matrix[cluster_matrix.medoid]
        edges = sorted(edges, key=lambda x: edges[x], reverse=True)
        cluster_matrix.reorder(edges)
        clusters[cluster] = cluster_matrix

    return sorted([c for c in clusters.values()], reverse=True)


def lloyds(matrix, medoids, eps=None, max_iterations=100):
    """Cluster the matrix using Lloyd's algorithm to the indicated
    threshold.

    :param matrix:
    :type matrix: SymMatrix
    :param medoids:
    :type medoids:
    :param eps:
    :type eps: float
    :param max_iterations:
    :type max_iterations: int
    :return: clusters
    """
    pass
