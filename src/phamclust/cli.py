"""Cluster phage genomes using gene content similarity-based metrics."""

import argparse
import pathlib
from multiprocessing import cpu_count

from phamclust.metrics import *


COLORS = "red,yellow,green"
CPUS = cpu_count()
EPILOG = """
Available metrics:

    Acronym     Name                                Reference
(1) gcs         gene content similarity             \
https://doi.org/10.1038/nmicrobiol.2017.112
(2) jc          jaccard coefficient                 \
https://doi.org/10.1111/j.1469-8137.1912.tb05611.x
(3) pocp        percentage of conserved proteins    \
https://doi.org/10.1128/JB.01688-14
(4) af          alignment fraction                  \
https://doi.org/10.1093/nar/gkv657
(5) aai         average aminoacid identity          \
https://doi.org/10.1073/pnas.0409727102
(6) peq         proteomic equivalence quotient      \
https://doi.org/10.1128/msystems.00443-23
"""
LINKAGES = {"single", "average", "complete"}
METRICS = {"gcs":  gene_content_similarity,
           "jc":   jaccard_coefficient,
           "pocp": percentage_of_conserved_proteins,
           "af":   alignment_fraction,
           "aai":  average_aminoacid_identity,
           "peq":  proteomic_equivalence_quotient}
K_MIN = 6                       # Minimum size of cluster to sub-cluster
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
                                description=__doc__, epilog=EPILOG,
                                formatter_class=argparse.RawTextHelpFormatter)

    p.add_argument("infile", type=pathlib.Path,
                   help=f"path to a TSV file mapping genomes to phams and "
                        f"translations")
    p.add_argument("outdir", type=pathlib.Path,
                   help=f"path to which output files should be written")
    p.add_argument("-g", "--genome-dir", action="store_true",
                   help=f"interpret `infile` as a directory of genome FASTA "
                        f"files instead of TSV")

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
    c.add_argument("-m", "--metric",
                   type=str, choices=METRICS, default=METRIC, metavar="",
                   help=f"relatedness index to use for pairwise genome "
                        f"comparisons [default: %(default)s]")

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
