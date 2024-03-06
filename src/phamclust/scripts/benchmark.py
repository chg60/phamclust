"""Script to benchmark CPU scaling of PhamClust distance matrix
calculations."""

import sys
import argparse
import datetime
import pathlib
import random

from phamclust.cli import METRICS
from phamclust.matrix import matrix_de_novo, matrix_to_squareform
from phamclust.scripts.phamclust import load_genomes_from_tsv

# Common values that might be found/used on modern commodity hardware
CPUS = [1, 2, 4, 6, 8, 12, 16]

SAMPLE_SIZES = [2000, 1000, 500]

# Default: run each test this many times to get representative outcomes
ITER = 3

# Random seed
SEED = 42


def parse_args():
    """Parse command line arguments."""
    p = argparse.ArgumentParser(description=__doc__, prog='phamclust-benchmark')

    p.add_argument("infile", type=pathlib.Path,
                   help="path to input TSV file mapping phages to phams and "
                        "translations")
    p.add_argument("outdir", type=pathlib.Path,
                   help="path where benchmark files can be written")

    p.add_argument("-m", "--metric",
                   type=str, default="peq", choices=METRICS.keys(), metavar="",
                   help=f"metric to use for pairwise distance calculations "
                        f"[default: %(default)s]")
    p.add_argument("-i", "--iterations",
                   type=int, default=ITER, metavar="",
                   help="number of times to run each test "
                        "[default: %(default)s]")
    p.add_argument("-s", "--seed",
                   type=int, default=SEED, metavar="",
                   help="random seed for reproducibility "
                        "[default: %(default)s]")

    return p.parse_args()


def dump_genomes_to_tsv(genomes, filepath):
    """Parse the 2- or 3-column input file into a dictionary of
    genomes.

    :param genomes: list of genomes
    :type genomes: list[phamclust.genome.Genome]
    :param filepath: the TSV file to dump genomes to
    :type filepath: str | pathlib.Path | os.PathLike[str]
    :return: filepath
    :rtype: pathlib.Path
    """
    with open(filepath, "w") as pham_writer:
        for genome in genomes:
            for pham, translations in genome:
                for translation in translations:
                    pham_writer.write(f"{genome.name}\t{pham}\t{translation}\n")

    return filepath


def main():
    """Commandline entrypoint."""
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    args = parse_args()

    # Get outdir and infile
    infile = args.infile
    outdir = args.outdir

    # Get metric
    metric = METRICS[args.metric]

    # Set random seed
    random.seed(args.seed)

    # Create outdir if not exists
    if not outdir.is_dir():
        outdir.mkdir(parents=True)

    # Load genomes from infile
    genomes = load_genomes_from_tsv(infile)

    # Sample different datasets from genomes
    for k in SAMPLE_SIZES:
        sample_genomes = random.sample(population=genomes, k=k)
        sample_file = outdir.joinpath(f"{k}_genomes.tsv")
        dump_genomes_to_tsv(sample_genomes, sample_file)

        for n in CPUS:
            for i in range(args.iterations):
                filename = outdir.joinpath(f"{k}_genomes-{n}_cpus-iter_{i}-"
                                           f"pairwise_{metric}_distances.tsv")
                t_start = datetime.datetime.now()

                matrix = matrix_de_novo(sample_genomes, metric, n)

                t_stop = datetime.datetime.now()
                matrix_to_squareform(matrix, filename)

                print(f"Genomes: {k}, CPUs: {n}, Iter: {i}, "
                      f"Elapsed: {str(t_stop - t_start)}")


if __name__ == "__main__":
    main()
