"""A program for assorting genomes into clusters and sub-clusters
based on shared gene content. Proteins evolve more slowly than the DNA
that encodes them, making this approach much more robust for comparing
the relatedness of bacteria between genera, or phages at almost any
genetic distance."""

import sys
import argparse
from pathlib import Path

from phamclust import clustering, similarity, multiprocess as mp
from phamclust.Genome import Genome
from phamclust.SymMatrix import matrix_de_novo

metrics = {"pocp": similarity.percentage_of_conserved_proteins,
           "gcs": similarity.gene_content_similarity,
           "jaccard": similarity.jaccard_similarity_coefficient}

CWD = Path().cwd()
METRIC = "pocp"
THRESH = 35


def parse_args():
    """Parse commandline arguments."""
    p = argparse.ArgumentParser(description=__doc__, prog="phamclust")

    p.add_argument("infile", type=Path,
                   help="path to tab-delimited file containing genome-pham "
                        "mappings")
    p.add_argument("outdir", type=Path, default=CWD,
                   help=f"directory where output files should be written "
                        f"[default: {CWD}]")

    # Clustering parameters
    c = p.add_argument_group("clustering parameters (optional)")
    c.add_argument("-m", "--metric", default=METRIC, type=str,
                   choices=metrics.keys(), metavar="",
                   help=f"metric to use for calculating pairwise similarity "
                        f"[default: {METRIC}]")
    c.add_argument("-t", "--threshold", default=THRESH, type=float, metavar="",
                   help=f"similarity threshold for a genome to join a cluster "
                        f"[default: {THRESH}%%]")
    c.add_argument("-s", "--subcluster", action="store_true",
                   help="break clusters into subclusters")

    p.add_argument("-c", "--cpus", default=mp.CPUS, type=int, metavar="",
                   help=f"number of CPU cores to use [default: {mp.CPUS}]")

    return p.parse_args()


def build_genomes(infile):
    """Parse genome-pham edges from the input file and populate a
    symmetric matrix with pairwise genome similarities.

    :param infile: tab-delimited file containing genome-pham edges
    :type infile: Path
    :return: matrix
    """
    genomes = dict()
    with open(infile, "r") as edge_reader:
        for row in edge_reader:
            name, pham = row.rstrip().split()
            if name not in genomes:
                genomes[name] = Genome(name)

            genomes[name].add(pham)

    return genomes


def get_subclusters(clusters):
    subclusters = [[]] * len(clusters)

    for i, cluster in enumerate(clusters):
        # If cluster too small to robustly determine eps, continue
        if len(cluster) < 5:
            continue

        # Otherwise - check eps to decide whether to subdivide
        eps = max([50, round(cluster.epsilon)])
        if eps >= 90:
            continue

        subclusters[i].extend(clustering.dbscan(cluster, eps))
        cluster_order = list()
        for j, subcluster in enumerate(subclusters[i]):
            cluster_order.extend(subcluster.node_names)
        cluster.reorder(cluster_order)

    return subclusters


def get_clusters(infile, threshold, func=metrics[METRIC], cpus=1):
    """Assemble clusters and subclusters.

    :param infile: tab-delimited file containing genome-pham edges
    :type infile: Path
    :param threshold: threshold to use for clustering
    :type threshold: int
    :param func: function to use for calculating similarity
    :type func: function
    :param cpus: number of CPU cores to use
    :type cpus: int
    :return: clusters
    """
    genomes = build_genomes(infile)

    matrix = matrix_de_novo(genomes, func, cpus)

    return clustering.dbscan(matrix, threshold)


def main():
    """Commandline interface for this module."""
    if len(sys.argv) == 1:
        sys.argv.append("-h")

    namespace = parse_args()

    infile = namespace.infile
    if not infile.is_file():
        print(f"input file ({infile}) is not a valid file - exiting...")
        sys.exit(1)

    outdir = namespace.outdir
    if not outdir.is_dir():
        print(f"output directory ({outdir}) does not exist - making it...")
        outdir.mkdir(parents=True)
    else:
        print(f"output directory ({outdir}) exists - overwriting contents...")
        response = input(f"press 'Enter' to continue, any other key to abort")
        if response:
            print(f"received abort signal - exiting...")
            sys.exit(1)

    metric = metrics[namespace.metric]
    threshold = namespace.threshold
    cpus = namespace.cpus
    do_sub = namespace.subcluster

    clusters = get_clusters(infile, threshold, metric, cpus)
    if do_sub:
        subclusters = get_subclusters(clusters)
    else:
        subclusters = [[]] * len(clusters)

    for i, group in enumerate(zip(clusters, subclusters)):
        _cluster, _subclusters = group

        _cluster.write_clustalo(outdir.joinpath(f"cluster_{i+1}.txt"))
        if not _subclusters or len(_subclusters) == 1:
            continue

        for j, _subcluster in enumerate(_subclusters):
            outfile = outdir.joinpath(f"cluster_{i+1}.{j+1}.txt")
            _subcluster.write_clustalo(outfile)


if __name__ == "__main__":
    main()
