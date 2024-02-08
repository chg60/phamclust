"""Commandline entrypoint to this module."""

import datetime
import hashlib
import logging
import pathlib
import shutil
import sys

from phamclust.cli import parse_args, METRICS
from phamclust.clustering import hierarchical
from phamclust.genome import Genome
from phamclust.heatmap import draw_heatmap
from phamclust.matrix import matrix_de_novo, matrix_to_squareform, \
    matrix_from_squareform

LOG_STR_FMT = "phamclust: %(asctime)s.%(msecs)03d: %(levelname)s: %(message)s"
LOG_TIME_FMT = "%H:%M:%S"


def load_genomes_from_tsv(filepath):
    """Parse the 2- or 3-column input file into a dictionary of
    genomes.

    :param filepath: the TSV file from which to parse genomes
    :type filepath: str | pathlib.Path | os.PathLike[str]
    :return: genomes
    :rtype: list[Genome]
    """
    genomes = dict()
    with open(filepath, "r") as pham_reader:
        for row in pham_reader:
            row = row.rstrip().split()

            if len(row) == 3:
                name, pham, translation = row
            elif len(row) == 2:
                (name, pham), translation = row, "M"
            else:
                raise ValueError(f"input file must either 2 or 3 columns")

            if name not in genomes:
                genomes[name] = Genome(name)

            genomes[name].add(pham, translation)

    return [genome for _, genome in genomes.items()]


def load_genomes_from_fasta_dir(filepath):
    """Load genomes from a directory of FASTA files. Each file must
    contain all the translations from a single genome, and the FASTA
    headers must be structured as:

    >name={name}|pham={pham}|n={n}

    where name is the genome identifier, pham is the pham identifier
    for the translation on the next line, and n is a counter for each
    paralogous gene copy in the genome in question.

    :param filepath: path to a directory containing genome FASTA(s)
    :type filepath: str | pathlib.Path | os.PathLike[str]
    :return: genomes
    :rtype: list[Genome]
    """
    genomes = dict()

    for f in filepath.iterdir():
        if f.suffix not in (".fasta", ".faa", ".fa"):
            logging.debug(f"{f} does not appear to be FASTA - skipping")
            continue
        elif f.is_dir():
            logging.debug(f"{f} appears to be a directory - skipping")
            continue

        name = f.stem
        if name not in genomes:
            genomes[name] = Genome(name)
            genomes[name].load(f)
        else:
            raise ValueError(f"duplicate genome name detected {name}")

    return [genome for _, genome in genomes.items()]


def _hash_genomes(genomes):
    """Create an MD5 hashsum from a list of genomes.

    This function assumes that genomes are sorted as desired. Changing
    the sort order WILL change the MD5 hashsum.

    Facilitates re-use of temporary files on subsequent runs of the
    same dataset.

    :param genomes: the genomes from which to compute a hashsum
    :type genomes: list[Genome]
    :return: md5_sum
    :rtype: str
    """
    hashsum = hashlib.new("md5")

    # Always traverse the genomes in ascending name order
    for genome in genomes:
        hashsum.update(str(genome).encode())

    return hashsum.hexdigest()


def check_matrix_integrity(matrix):
    """Verify that a matrix has been filled correctly.

    Given a matrix size `n`, verifies that the number of edges is
    correct. Other checks include verifying that the edges are not missing
    weights and that the diagonal sum matches the expectation for
    distance/identity matrices (whichever the matrix claims to be).

    :param matrix: the matrix whose integrity is to be verified
    :type matrix: phamclust.matrix.SymMatrix
    :return: is_ok
    """
    cell_exp = int(len(matrix) * (len(matrix) - 1) / 2)
    diag_exp = (1.0 - matrix.is_distance) * len(matrix)
    cell_cnt, empty_cnt, diag_sum = 0, 0, 0
    for source, target, weight in matrix:
        if source == target:
            diag_sum += weight
        else:
            cell_cnt += 1

        if weight is None:
            empty_cnt += 1

    cell_diff = cell_cnt - cell_exp
    diag_diff = diag_sum - diag_exp
    none_diff = empty_cnt - 0

    return cell_diff, diag_diff, none_diff


def phamclust(infile, outdir, is_genome_dir, metric, nr_distance, nr_linkage,
              clu_distance, clu_linkage, sub_distance, sub_linkage, k_min,
              no_sub, cpus, rm_tmp, debug):
    """Main program for phamclust.

    :param infile: input TSV file mapping phage-to-pham-to-translation
    :type infile: pathlib.Path
    :param outdir: output directory where files should be written
    :type outdir: pathlib.Path
    :param is_genome_dir: indicate that infile should instead be
        interpreted as a directory of genome FASTA files
    :type is_genome_dir: bool
    :param metric: function to use for distance calculations
    :type metric: str
    :param nr_distance: distance threshold for clustering very similar
        genomes
    :type nr_distance: float
    :param nr_linkage: hierarchical clustering linkage type to use for
        clustering very similar genomes
    :type nr_linkage: str
    :param clu_distance: distance threshold for clustering less similar
        genomes
    :type clu_distance: float
    :param clu_linkage: hierarchical clustering linkage type to use for
        clustering less similar genomes
    :type clu_linkage: str
    :param sub_distance: distance threshold for sub-clustering genome
        clusters
    :type sub_distance: float
    :param sub_linkage: hierarchical clustering linkage type to use for
        sub-clustering genome clusters
    :type sub_linkage: str
    :param k_min: minimum size of cluster to perform sub-clustering on
    :type k_min: int
    :param no_sub: indicate that sub-clustering should not be performed
    :type no_sub: bool
    :param cpus: number of CPUs to use for distance calculations
    :type cpus: int
    :param rm_tmp: remove temporary files when done
    :type rm_tmp: bool
    :param debug: indicate if logging should be set to DEBUG level
    :type debug: bool
    """
    # Sanity check - if nr_dist not lower than clu_dist, use 0.0
    if nr_distance >= clu_distance:
        nr_distance = 0.0

    logging.info(f"=======================")
    logging.info(f" 0: runtime parameters ")
    logging.info(f"=======================")
    logging.info(f"infile:     {infile}")
    logging.info(f"outdir:     {outdir}")
    logging.info(f"debug:      {debug}")
    logging.info(f"subcluster: {not no_sub}")
    logging.info(f"remove tmp: {rm_tmp}")
    logging.info(f"sub dist:   {sub_distance}")
    logging.info(f"sub link:   {sub_linkage}")
    logging.info(f"clu dist:   {clu_distance}")
    logging.info(f"clu link:   {clu_linkage}")
    logging.info(f"nr dist:    {nr_distance}")
    logging.info(f"nr link:    {nr_linkage}")
    logging.info(f"metric:     {metric}")
    logging.info(f"cpus:       {cpus}")

    # Parse genomes from TSV-formatted infile
    logging.info(f"====================")
    logging.info(f" 1: parsing genomes ")
    logging.info(f"====================")
    if is_genome_dir:
        genomes = load_genomes_from_fasta_dir(infile)
        logging.info(f"loaded {len(genomes)} genomes from input directory")
    else:
        genomes = load_genomes_from_tsv(infile)
        logging.info(f"loaded {len(genomes)} genomes from input TSV")

    genomes.sort(key=lambda x: x.name)
    hashsum = _hash_genomes(genomes)
    logging.info(f"md5 hashsum is {hashsum}")

    # Create temporary directory (or re-use) from genome hashsum
    tmpdir = outdir.joinpath(f"{hashsum}.tmp")
    if not tmpdir.is_dir():
        tmpdir.mkdir()
    logging.info(f"using temp directory {tmpdir.name}")

    # Write genomes to temporary directory as FASTA
    tmp_genomes = tmpdir.joinpath(f"01_genomes")
    if not tmp_genomes.is_dir():
        tmp_genomes.mkdir()

    for genome in genomes:
        genome_fasta = tmp_genomes.joinpath(f"{genome.name}.fasta")
        if not genome_fasta.is_file():
            genome.save(genome_fasta)

    logging.info(f"genomes stashed in {tmpdir.name}/{tmp_genomes.name}")

    # Always recycle pairwise matrices if available
    logging.info(f"==========================")
    logging.info(f" 2: build distance matrix ")
    logging.info(f"==========================")
    logging.info(f"selected metric: {metric}")

    # Write distance matrix to temporary directory as FASTA
    tmp_distmats = tmpdir.joinpath(f"02_distmats")
    if not tmp_distmats.is_dir():
        tmp_distmats.mkdir()

    dist_mat_file = tmp_distmats.joinpath(f"{metric}_distance_matrix.tsv")
    if not dist_mat_file.is_file():
        logging.info(f"cached distance matrix not found - computing de novo")
        _start = datetime.datetime.now()
        dist_mat = matrix_de_novo(genomes, METRICS[metric], cpus)
        _elapsed = datetime.datetime.now() - _start
        logging.info(f"computed distance matrix in {str(_elapsed)}")
        logging.info(f"caching distance matrix so it can be re-used")
        matrix_to_squareform(dist_mat, dist_mat_file, lower_triangle=True)
    else:
        logging.info(f"found cached distance matrix - importing it")
        _start = datetime.datetime.now()
        dist_mat = matrix_from_squareform(dist_mat_file)
        _elapsed = datetime.datetime.now() - _start
        logging.info(f"loaded distance matrix in {str(_elapsed)}")

    # Make sure the matrix is a distance matrix
    if not dist_mat.is_distance:
        logging.info(f"matrix is not a distance matrix - flipping it")
        dist_mat.invert()

    # Check matrix integrity
    status = check_matrix_integrity(dist_mat)
    if not any(status):
        logging.info(f"matrix passed all integrity checks")
    else:
        logging.error(f"matrix failed the following integrity check(s):")
        if status[0] > 0:
            logging.error(f"found {status[0]} too many edges")
        elif status[0] < 0:
            logging.error(f"found {0 - status[0]} too few edges")
        if status[1] != 0:
            logging.error(f"diagonal sum {status[0]} is wrong")
        if status[2] != 0:
            logging.error(f"found {status[2]} unfilled edges")

        print(f"matrix validation failed - check log for details")
        sys.exit(1)

    # If we got here, matrix is OK - begin clustering
    logging.info(f"====================")
    logging.info(f" 3: cluster genomes ")
    logging.info(f"====================")
    logging.info(f"grouping highly redundant genomes with distance <= "
                 f"{nr_distance} by {nr_linkage} linkage")
    # Seed clusters with groups of phages that have low levels of divergence
    clu_mats = hierarchical(dist_mat, eps=nr_distance, linkage=nr_linkage)

    # Map the seed cluster representatives to their clusters
    seed_map = {x.medoid[0]: x for x in clu_mats}

    # Merge seed clusters with high aligned fractions
    repr_mat = dist_mat.extract_submatrix(list(seed_map.keys()))
    logging.info(f"found {len(repr_mat)} groups of similar genomes")
    logging.info(f"clustering non-redundant genomes with distance <= "
                 f"{clu_distance} by {clu_linkage} linkage")
    clu_mats = hierarchical(repr_mat, eps=clu_distance, linkage=clu_linkage)
    for i, clu_mat in enumerate(clu_mats):
        all_nodes = list()
        for repr_node in clu_mat.nodes:
            all_nodes.extend(seed_map[repr_node].nodes)
        clu_mats[i] = dist_mat.extract_submatrix(all_nodes)

    single_mats = [clu_mat for clu_mat in clu_mats if len(clu_mat) == 1]
    clu_mats = [clu_mat for clu_mat in clu_mats if len(clu_mat) > 1]
    logging.info(f"found {len(clu_mats)} clusters and {len(single_mats)} "
                 f"singletons")

    # Write genomes to temporary directory as FASTA
    tmp_clusters = tmpdir.joinpath(f"03_clusters")
    if not tmp_clusters.is_dir():
        tmp_clusters.mkdir()

    # Begin sub-clustering
    logging.info(f"========================")
    logging.info(f" 4: sub-cluster genomes ")
    logging.info(f"========================")
    for i, clu_mat in enumerate(sorted(clu_mats, reverse=True)):
        logging.info(f"cluster {i + 1} has {len(clu_mat)} nodes")

        # Make a toplevel directory for this cluster
        cluster_dir = tmp_clusters.joinpath(f"cluster_{i+1}")
        if cluster_dir.is_dir():
            logging.debug(f"removing old cluster directory {cluster_dir}")
            shutil.rmtree(cluster_dir)

        logging.debug(f"creating new cluster directory {cluster_dir}")
        cluster_dir.mkdir()

        cluster_nodes = set(clu_mat.nodes)
        cluster_genomes = [x for x in genomes if x.name in cluster_nodes]
        cluster_matfile = cluster_dir.joinpath(f"{metric}_similarity.tsv")
        cluster_heatmap = cluster_dir.joinpath(f"{metric}_heatmap.pdf")
        cluster_heatmap2 = cluster_dir.joinpath(f"{metric}_heatmap.html")

        # Make a genome directory for this cluster, move genome files here
        genome_dir = cluster_dir.joinpath(f"genomes")
        genome_dir.mkdir()
        for genome in cluster_genomes:
            genome_file = genome_dir.joinpath(f"{genome.name}.faa")
            genome.save(genome_file)
            logging.debug(f"wrote genome to {genome_file}")

        if no_sub or len(clu_mat) < k_min:
            logging.debug(f"not sub-clustering {len(clu_mat)} genomes")

            # Re-order matrix in order of proximity to the medoid; convert to
            # similarity matrix
            order = [clu_mat.medoid[0]]
            order.extend(clu_mat.nearest_neighbors(order[0], threshold=1.0))
            clu_mat.reorder(order)
            clu_mat.invert()
            matrix_to_squareform(clu_mat, cluster_matfile)
            draw_heatmap(clu_mat, midpoint=0.5, filename=cluster_heatmap)
            draw_heatmap(clu_mat, midpoint=0.5, filename=cluster_heatmap2)
            continue
        else:
            sub_mats = hierarchical(clu_mat, eps=sub_distance,
                                    linkage=sub_linkage)
            sub_mats = sorted(sub_mats, reverse=True)
            order = list()
            for j, sub_mat in enumerate(sub_mats):
                sub_matfile = cluster_dir.joinpath(
                    f"subcluster_{j+1}_similarity.tsv")
                _ord = [sub_mat.medoid[0]]
                _ord.extend(sub_mat.nearest_neighbors(_ord[0], threshold=1.0))
                order.extend(_ord)
                sub_mat.reorder(_ord)
                sub_mat.invert()
                matrix_to_squareform(sub_mat, sub_matfile)
            clu_mat.reorder(order)
            clu_mat.invert()
            matrix_to_squareform(clu_mat, cluster_matfile)
            draw_heatmap(clu_mat, midpoint=0.5, filename=cluster_heatmap)
            draw_heatmap(clu_mat, midpoint=0.5, filename=cluster_heatmap2)
            continue

    # Make a toplevel directory for singletons, with nested genome dir
    if single_mats:
        single_dir = tmp_clusters.joinpath(f"singletons")
        if single_dir.is_dir():
            logging.debug(f"removing old singleton directory {single_dir}")
            shutil.rmtree(single_dir)

        logging.debug(f"creating new singleton directory {single_dir}")
        single_dir.mkdir()
        genome_dir = single_dir.joinpath("genomes")
        genome_dir.mkdir()

        for i, single_mat in enumerate(single_mats):
            node = single_mat.nodes[0]
            genome = [x for x in genomes if x.name == node][0]

            genome_file = genome_dir.joinpath(f"{genome.name}.faa")
            genome.save(genome_file)
            logging.debug(f"wrote genome to {genome_file}")

    # Write full dataset heatmap
    logging.info(f"========================")
    logging.info(f" 5: draw dataset heatmap")
    logging.info(f"========================")
    logging.info(f"putting matrix in cluster order")
    dataset_order = list()
    for clu_mat in sorted(clu_mats, reverse=True):
        dataset_order.extend(clu_mat.nodes)
    for single_mat in single_mats:
        dataset_order.extend(single_mat.nodes)

    dataset_heatmap = tmp_clusters.joinpath(f"{metric}_heatmap.html")
    dist_mat.reorder(dataset_order)

    logging.info(f"convert to similarity matrix for easier visualization")
    dist_mat.invert()

    logging.info(f"drawing heatmap and saving to {dataset_heatmap.name}")
    draw_heatmap(dist_mat, midpoint=0.5, filename=dataset_heatmap)

    # Copy output files to output directory
    logging.info(f"========================")
    logging.info(f" 6: move output files   ")
    logging.info(f"========================")
    # if outdir.is_dir():
    #     logging.info(f"removing existing output directory {outdir}")
    #     shutil.rmtree(outdir)

    logging.info(f"moving output files from temporary directory to {outdir}")
    shutil.copytree(tmp_clusters, outdir, dirs_exist_ok=True)

    # Cleanup temporary directory if asked
    if rm_tmp:
        logging.info(f"cleaning up temporary files in {tmpdir}")
        shutil.rmtree(tmpdir)


def main():
    """Commandline entry point."""
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    args = parse_args()

    # Check infile and outdir
    if args.genome_dir and not args.infile.is_dir():
        print(f"genome directory '{args.infile}' does not exist")
        sys.exit(1)
    elif not args.genome_dir and not args.infile.is_file():
        print(f"input TSV '{args.infile}' does not exist")
        sys.exit(1)

    if not args.outdir.is_dir():
        args.outdir.mkdir(parents=True)

    # Set logging context and log runtime parameters
    logfile = args.outdir.joinpath("phamclust.log")
    log_level = logging.INFO
    if args.debug:
        log_level = logging.DEBUG

    logging.basicConfig(filename=logfile, filemode="w", level=log_level,
                        format=LOG_STR_FMT, datefmt=LOG_TIME_FMT)
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    # Run PhamClust!
    phamclust(infile=args.infile,
              outdir=args.outdir,
              is_genome_dir=args.genome_dir,
              metric=args.metric,
              nr_distance=round(1.0 - args.nr_thresh, 6),
              nr_linkage=args.nr_linkage,
              clu_distance=round(1.0 - args.clu_thresh, 6),
              clu_linkage=args.clu_linkage,
              sub_distance=round(1.0 - args.sub_thresh, 6),
              sub_linkage=args.sub_linkage,
              k_min=max([1, args.k_min]),
              cpus=args.threads,
              no_sub=args.no_sub,
              rm_tmp=args.remove_tmp,
              debug=args.debug)


if __name__ == "__main__":
    main()
