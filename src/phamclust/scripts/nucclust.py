"""PhamClust pipeline, but using BLASTN to calculate inter-genomic distances."""

import argparse
import csv
import datetime
import logging
import pathlib
import shlex
import shutil
import sys
from multiprocessing import cpu_count
from subprocess import Popen, PIPE

import joblib

from phamclust.clustering import hierarchical_clustering
from phamclust.heatmap import draw_heatmap, CSS_COLORS
from phamclust.matrix import matrix_to_adjacency, \
    matrix_from_squareform, matrix_to_squareform, SymMatrix

# BLAST fields to use
BLASTN_FIELDS = ['qseqid', 'sseqid', 'nident', 'length', 'qstart', 'qend',
                 'sstart', 'send', 'qlen', 'slen', 'evalue', 'bitscore']

# Heatmap settings
COLORS = "red,yellow,green"
CPUS = cpu_count()

# Logging settings
LOG_STR_FMT = "nucclust: %(asctime)s.%(msecs)03d: %(levelname)s: %(message)s"
LOG_TIME_FMT = "%H:%M:%S"

# Clustering settings
K_MIN = 6                       # Minimum size of cluster to sub-cluster
LINKAGES = {"single", "average", "complete"}    # Available linkage types
METRIC = "peq"                  # Default metric
NR_THRESH = 0.75                # Default 1st iteration cluster threshold
NR_LINKAGE = "complete"         # Default 1st iteration cluster linkage
CLU_THRESH = 0.25               # Default 2nd iteration cluster threshold
CLU_LINKAGE = "average"         # Default 2nd iteration cluster linkage
SUB_THRESH = 0.6                # Default subclustering threshold
SUB_LINKAGE = "single"          # Default subclustering linkage


def parse_args():
    """Parse commandline arguments."""
    p = argparse.ArgumentParser(prog="phamclust",
                                description=__doc__,
                                formatter_class=argparse.RawTextHelpFormatter)

    p.add_argument("indir", type=pathlib.Path,
                   help=f"path to a directory containing nucleotide FASTA "
                        f"files")
    p.add_argument("outdir", type=pathlib.Path,
                   help=f"path to which output files should be written")

    c = p.add_argument_group("clustering arguments:")
    c.add_argument("-k", "--k-min",
                   type=int, default=K_MIN, metavar="",
                   help=f"minimum cluster size to perform subclustering "
                        f"[default: %(default)s]")
    c.add_argument("-s", "--sub-thresh",
                   type=float, default=SUB_THRESH, metavar="",
                   help=f"similarity threshold to use for sub-clustering "
                        f"[default: %(default)s]")
    c.add_argument("-sl", "--sub-linkage",
                   type=str, choices=LINKAGES, default=SUB_LINKAGE, metavar="",
                   help=f"linkage type to use for sub-clustering "
                        f"[default: %(default)s]")
    c.add_argument("-c", "--clu-thresh",
                   type=float, default=CLU_THRESH, metavar="",
                   help=f"similarity threshold to use for clustering "
                        f"[default: %(default)s]")
    c.add_argument("-cl", "--clu-linkage",
                   type=str, choices=LINKAGES, default=CLU_LINKAGE, metavar="",
                   help=f"linkage type to use for clustering "
                        f"[default: %(default)s]")
    c.add_argument("-nr", "--nr-thresh",
                   type=float, default=NR_THRESH, metavar="",
                   help=f"similarity threshold above which to pre-group "
                        f"very similar genomes that must be clustered together "
                        f"[default: %(default)s]")
    c.add_argument("-nl", "--nr-linkage",
                   type=str, choices=LINKAGES, default=NR_LINKAGE, metavar="",
                   help=f"linkage type to use for pre-grouping very similar "
                        f"genomes [default: %(default)s]")

    h = p.add_argument_group("heatmap arguments:")
    h.add_argument("-hc", "--heatmap-colors",
                   type=str, default=COLORS, metavar="",
                   help=f"comma-separated list of 2 or 3 colors to use in "
                        f"heatmaps [default: %(default)s]")
    h.add_argument("-hm", "--heatmap-midpoint",
                   type=float, default=0.5, metavar="",
                   help=f"midpoint to use for color gradient in heatmaps "
                        f"[default: %(default)s]")

    p.add_argument("-d", "--debug", action="store_true",
                   help=f"increase verbosity of logging for debug purposes")
    p.add_argument("-n", "--no-sub", action="store_true",
                   help="do not perform sub-clustering")
    p.add_argument("-r", "--remove-tmp", action="store_true",
                   help=f"remove temporary files (not recommended if repeated "
                        f"runs are planned on the same dataset)")
    p.add_argument("-t", "--threads",
                   type=int, default=CPUS, metavar="",
                   help=f"number of CPU cores to use [default: %(default)s]")

    return p.parse_args()


class BlastError(Exception):
    """An error raised due to a problem running BLAST."""
    pass


def blastn(query, subject):
    """Run blastn on a pair of nucleotide FASTA files and return the stdout.

    :param query: the query FASTA file
    :type query: pathlib.Path
    :param subject: the subject FASTA file
    :type subject: pathlib.Path
    :return: the stdout from blastn
    """
    fields = " ".join(BLASTN_FIELDS)
    command = f"blastn -query {query} -subject {subject} -outfmt '6 " \
              f"{fields}' -task blastn -evalue 0.001 -dust no -culling_limit 1"

    with Popen(shlex.split(command), stdout=PIPE, stderr=PIPE) as proc:
        stdout = proc.stdout.read().decode("utf-8")
        stderr = proc.stdout.read().decode("utf-8")

    if stderr:
        raise BlastError(stderr)

    return stdout


def blastn_multiple(query, subjects, d):
    """Run blastn between a single query and multiple subjects and
    write the results to a file in `d` named after the query.

    :param query: the query FASTA file
    :type query: pathlib.Path
    :param subjects: a list of subject FASTA files
    :type subjects: list[pathlib.Path]
    :param d: the output directory to write results to
    :type d: pathlib.Path
    :return: the output filepath
    """
    outname = d.joinpath(f"{query.stem}.tsv")

    for subject in subjects:
        with open(outname, "a") as blast_handle:
            blast_handle.write(blastn(query, subject))

    return outname


class Genome:
    """A genome object representing a single nucleotide sequence."""
    def __init__(self, filepath):
        self.filepath = filepath

    @property
    def name(self):
        """The name of the genome."""
        return self.filepath.stem


def matrix_de_novo(genomes, cpus, tmp, as_distance=True):
    """Compute a pairwise matrix of distances between genomes.

    :param genomes: a list of genomes to compare
    :type genomes: list[Genome]
    :param cpus: the number of CPUs to use for distance calculations
    :type cpus: int
    :param tmp: the temporary directory to write intermediate files to
    :type tmp: pathlib.Path
    :param as_distance: indicate if the matrix should be a distance matrix
    :type as_distance: bool
    :return: a pairwise distance matrix
    """
    # >>> Begin sanity checks >>>
    if len(genomes) == 0:
        raise ValueError(f"need at least 1 genome to construct matrix de novo")

    # Calculate the number of unique genome pairs
    n_pairs = int((len(genomes) * (len(genomes) - 1) / 2) + len(genomes))
    logging.debug(f"{len(genomes)} genomes -> {n_pairs} edges including "
                  f"diagonal")

    # If fewer genome pairs than CPUs, set CPU count to number of genome pairs
    if n_pairs < cpus:
        logging.info(f"small dataset - reducing # CPUs to {n_pairs}")
        cpus = n_pairs
    # <<< End sanity checks <<<

    # <<< Initialize distance matrix >>>
    matrix = SymMatrix(nodes=[g.name for g in genomes], is_distance=as_distance)

    # Initialize parallel runner - disable memory mapping
    runner = joblib.Parallel(n_jobs=cpus, return_as="generator_unordered",
                             max_nbytes=None)

    # Set up and dispatch the parallel jobs
    batch, temp_outs = list(), list()
    for query in genomes:
        outfile = tmp.joinpath(f"{query.name}.tsv")
        if outfile.is_file():
            temp_outs.append(outfile)
            continue
        batch.append((query.filepath, [t.filepath for t in genomes], tmp))
    temp_outs.extend(runner(joblib.delayed(blastn_multiple)(*b) for b in batch))

    # Read BLASTN data into a dictionary mapping all hit data between observed
    # source-target pairs
    # BLASTN_FIELDS = ['qseqid', 'sseqid', 'nident', 'length', 'qstart', 'qend',
    #                  'sstart', 'send', 'qlen', 'slen', 'evalue' 'bitscore']
    blastn_map = dict()
    for temp_out in temp_outs:
        with open(temp_out, "r") as temp_reader:
            reader = csv.DictReader(temp_reader, delimiter="\t",
                                    fieldnames=BLASTN_FIELDS)
            # Iterate over rows (no header, so no need to skip first row)
            for row in reader:
                # source, target, *(data) = row.rstrip().split("\t")
                source, target = row.pop("qseqid"), row.pop("sseqid")
                for key, value in row.items():
                    try:
                        row[key] = int(value)
                    except ValueError:
                        row[key] = float(value)

                if source in blastn_map:
                    if target in blastn_map[source]:
                        blastn_map[source][target].append(row)
                    else:
                        blastn_map[source][target] = [row]
                else:
                    blastn_map[source] = {target: [row]}

    # blastn_map is now a dict() mapping source genome names to target
    # genome names to a list of BLASTN data key-value pairs
    nodes = matrix.nodes
    for i, source in enumerate(nodes):
        for target in nodes[i:]:
            # source_hits result from blastn(source, target)
            source_hits = blastn_map[source].get(target, dict())

            # target_hits result from blastn(target, source)
            target_hits = blastn_map[target].get(source, dict())

            # We prefer the case where source_hits and target_hits are both
            # non-empty, but we can handle the case where one is empty
            if source_hits and target_hits:
                numerator = sum([hsp["nident"] for hsp in source_hits])
                numerator += sum([hsp["nident"] for hsp in target_hits])
                denominator = source_hits[0]["qlen"] + target_hits[0]["qlen"]
                weight = min([numerator / denominator, 1.0])
            # If we have source_hits but no target_hits, we can calculate the
            # weight as the sum of the nident values divided by the qlen of the
            # source genome
            elif source_hits:
                numerator = sum([hsp["nident"] for hsp in source_hits])
                denominator = source_hits[0]["qlen"]
                weight = min([numerator / denominator, 1.0])
            # If we have target_hits but no source_hits, we can calculate the
            # weight as the sum of the nident values divided by the qlen of the
            # target genome
            elif target_hits:
                numerator = sum([hsp["nident"] for hsp in target_hits])
                denominator = target_hits[0]["qlen"]
                weight = min([numerator / denominator, 1.0])
            # If we have no recorded hits between source and target, weight is 0
            else:
                weight = 0.0

            if as_distance:
                weight = 1.0 - weight

            matrix.set_weight(source, target, weight)

    return matrix


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



def nucclust(indir, outdir, nr_distance, nr_linkage, clu_distance,
             clu_linkage, sub_distance, sub_linkage, k_min, no_sub, colors,
             midpoint, cpus, rm_tmp, debug):
    """Main program for nucclust.

    :param indir: input directory containing nucleotide FASTA files
    :type indir: pathlib.Path
    :param outdir: output directory where files should be written
    :type outdir: pathlib.Path
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
    :param colors: colors to use in heatmaps
    :type colors: list[str]
    :param midpoint: midpoint of the color gradient in heatmaps
    :type midpoint: float
    :param cpus: number of CPUs to use for distance calculations
    :type cpus: int
    :param rm_tmp: remove temporary files when done
    :type rm_tmp: bool
    :param debug: indicate if logging should be set to DEBUG level
    :type debug: bool
    """
    metric = "nucleotide_similarity"

    # Sanity check - if nr_dist not lower than clu_dist, use 0.0
    if nr_distance >= clu_distance:
        nr_distance = 0.0

    logging.info(f"=======================")
    logging.info(f" 0: runtime parameters ")
    logging.info(f"=======================")
    logging.info(f"indir:      {indir}")
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
    logging.info(f"colors:     {','.join(colors)}")
    logging.info(f"midpoint:   {midpoint}")

    # Parse genomes from TSV-formatted infile
    logging.info(f"====================")
    logging.info(f" 1: parsing genomes ")
    logging.info(f"====================")
    genomes = [Genome(x) for x in indir.iterdir() if x.is_file() and x.suffix
               in {".fna", ".fa", ".fasta"}]
    logging.info(f"loaded {len(genomes)} genomes from input directory")

    genomes.sort(key=lambda x: x.name)

    # Create temporary directory
    tmpdir = outdir.joinpath(f"tempfiles.tmp")
    if not tmpdir.is_dir():
        tmpdir.mkdir()
    logging.info(f"using temp directory {tmpdir.name}")

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
        dist_mat = matrix_de_novo(genomes, cpus, tmp_distmats)
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
            logging.error(f"diagonal sum {status[1]} is wrong")
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
    clu_mats = hierarchical_clustering(dist_mat, eps=nr_distance,
                                       linkage=nr_linkage)

    # Map the seed cluster representatives to their clusters
    seed_map = {x.medoid[0]: x for x in clu_mats}

    # Merge seed clusters with high aligned fractions
    repr_mat = dist_mat.extract_submatrix(list(seed_map.keys()))
    logging.info(f"found {len(repr_mat)} groups of similar genomes")
    logging.info(f"clustering non-redundant genomes with distance <= "
                 f"{clu_distance} by {clu_linkage} linkage")
    clu_mats = hierarchical_clustering(repr_mat, eps=clu_distance,
                                       linkage=clu_linkage)
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
        cluster_dir = tmp_clusters.joinpath(f"cluster_{i + 1}")
        if cluster_dir.is_dir():
            logging.debug(f"removing old cluster directory {cluster_dir}")
            shutil.rmtree(cluster_dir)

        logging.debug(f"creating new cluster directory {cluster_dir}")
        cluster_dir.mkdir()

        cluster_nodes = set(clu_mat.nodes)
        cluster_genomes = [x for x in genomes if x.name in cluster_nodes]
        cluster_matfile = cluster_dir.joinpath(f"{metric}_similarity.tsv")
        cluster_heatmap = cluster_dir.joinpath(f"{metric}_heatmap.svg")
        cluster_heatmap2 = cluster_dir.joinpath(f"{metric}_heatmap.html")

        # Make a genome directory for this cluster, move genome files here
        genome_dir = cluster_dir.joinpath(f"genomes")
        genome_dir.mkdir()
        for genome in cluster_genomes:
            genome_file = genome_dir.joinpath(f"{genome.name}.faa")
            with open(genome_file, "w") as fasta_writer:
                with open(genome.filepath, "r") as fasta_reader:
                    [fasta_writer.write(line) for line in fasta_reader]

        if no_sub or len(clu_mat) < k_min:
            logging.debug(f"not sub-clustering {len(clu_mat)} genomes")

            # Re-order matrix in distance tree order
            clu_mat.reorder()
            clu_mat.invert()
            matrix_to_squareform(clu_mat, cluster_matfile)
            draw_heatmap(clu_mat, colors=colors, midpoint=midpoint,
                         filename=cluster_heatmap)
            draw_heatmap(clu_mat, colors=colors, midpoint=midpoint,
                         filename=cluster_heatmap2)
            continue
        else:
            sub_mats = hierarchical_clustering(clu_mat, eps=sub_distance,
                                               linkage=sub_linkage)
            sub_mats = sorted(sub_mats, reverse=True)
            order = list()
            for j, sub_mat in enumerate(sub_mats):
                sub_matfile = cluster_dir.joinpath(
                    f"subcluster_{j + 1}_similarity.tsv")
                # Can't re-order a matrix with only a single node, and order
                # doesn't matter for just 2 nodes
                if len(sub_mat) > 2:
                    sub_mat.reorder()
                order.extend(sub_mat.nodes)
                sub_mat.invert()
                matrix_to_squareform(sub_mat, sub_matfile)
            clu_mat.reorder(order)
            clu_mat.invert()
            matrix_to_squareform(clu_mat, cluster_matfile)
            draw_heatmap(clu_mat, colors=colors, midpoint=midpoint,
                         filename=cluster_heatmap)
            draw_heatmap(clu_mat, colors=colors, midpoint=midpoint,
                         filename=cluster_heatmap2)
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
            with open(genome_file, "w") as fasta_writer:
                with open(genome.filepath, "r") as fasta_reader:
                    [fasta_writer.write(line) for line in fasta_reader]

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

    dataset_html = tmp_clusters.joinpath(f"{metric}_heatmap.html")
    dataset_svg = tmp_clusters.joinpath(f"{metric}_heatmap.svg")
    dist_mat.reorder(dataset_order)

    logging.info(f"cast to similarity matrix for easier visualization")
    dist_mat.invert()

    if len(dist_mat) > 1500:
        logging.info(f"full pairwise matrix is too large to visualize")
    else:
        logging.info(f"drawing heatmap HTML and saving to {dataset_html}")
        draw_heatmap(dist_mat, colors=colors, midpoint=1.0 - clu_distance,
                     filename=dataset_html)
        logging.info(f"drawing heatmap SVG and saving to {dataset_svg}")
        draw_heatmap(dist_mat, colors=colors, midpoint=1.0 - clu_distance,
                     filename=dataset_svg)

    # Copy output files to output directory
    logging.info(f"========================")
    logging.info(f" 6: move output files   ")
    logging.info(f"========================")
    logging.info(f"removing contents from existing output directory {outdir}")
    for fp in outdir.iterdir():
        if fp.name == tmpdir.name or fp.suffix == ".log":
            continue
        if fp.is_file():
            fp.unlink()
        elif fp.is_dir():
            shutil.rmtree(fp)
        else:
            logging.warning(f"skip removal of unknown filetype {fp}")

    matrix_out = outdir.joinpath(f"pairwise_{metric}_similarities.tsv")
    logging.info(f"writing pairwise {metric} similarities to {matrix_out}")
    matrix_to_squareform(dist_mat, matrix_out)

    matrix_out = outdir.joinpath(f"pairwise_{metric}_adjacency.tsv")
    logging.info(f"writing pairwise {metric} adjacency to {matrix_out}")
    matrix_to_adjacency(dist_mat, matrix_out, skip_zero=True)

    logging.info(f"moving output files from temporary directory to {outdir}")
    shutil.copytree(tmp_clusters, outdir, dirs_exist_ok=True)
    shutil.rmtree(tmp_clusters)

    # Cleanup temporary directory if asked
    if rm_tmp:
        logging.info(f"cleaning up temporary files in {tmpdir}")
        shutil.rmtree(tmpdir)


def main():
    """Commandline entrypoint."""
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    args = parse_args()

    # Check infile and outdir
    if not args.indir.is_dir():
        print(f"genome directory '{args.indir}' does not exist or is not a "
              f"directory")
        sys.exit(1)

    if not args.outdir.is_dir():
        args.outdir.mkdir(parents=True)

    # Set logging context and log runtime parameters
    logfile = args.outdir.joinpath("nucclust.log")
    log_level = logging.INFO
    if args.debug:
        log_level = logging.DEBUG

    logging.basicConfig(filename=logfile, filemode="w", level=log_level,
                        format=LOG_STR_FMT, datefmt=LOG_TIME_FMT)
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    # sanity-check colors - must be a list of either 2 or 3 valid colors
    colors = args.heatmap_colors.split(",")
    if not 2 <= len(colors) <= 3:
        logging.error(f"expected either two or three colors, got {len(colors)}")
        sys.exit(1)
    if len(colors) == 2:
        logging.warning("2-color scale ignores `--heatmap-midpoint`")

    valid_colors = []
    for color in colors:
        if color not in CSS_COLORS:
            logging.error(f"unknown color specified: '{color}'")
            continue
        valid_colors.append(color)
    if len(valid_colors) != len(colors):
        logging.error(f"got {len(colors) - len(valid_colors)} unrecognized "
                      f"colors")
        logging.error(f"valid colors:")
        logging.error(f"{' '.join(CSS_COLORS)}")
        sys.exit(1)

    # Run NucClust!
    nucclust(indir=args.indir,
             outdir=args.outdir,
             nr_distance=round(1.0 - args.nr_thresh, 6),
             nr_linkage=args.nr_linkage,
             clu_distance=round(1.0 - args.clu_thresh, 6),
             clu_linkage=args.clu_linkage,
             sub_distance=round(1.0 - args.sub_thresh, 6),
             sub_linkage=args.sub_linkage,
             k_min=max([1, args.k_min]),
             colors=valid_colors,
             midpoint=round(args.heatmap_midpoint, 6),
             cpus=args.threads,
             no_sub=args.no_sub,
             rm_tmp=args.remove_tmp,
             debug=args.debug)


if __name__ == "__main__":
    main()
