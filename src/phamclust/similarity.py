"""Functions for calculating various similarity metrics between pairs
of genomes."""

from phamclust.Genome import Genome
from phamclust.statistics import average


def percentage_of_conserved_proteins(source, target):
    """Calculate the percentage of conserved proteins between the
    given source and target genomes.

    Qin et al. A Proposed Genus Boundary for the Prokaryotes Based on
    Genomic Insights.
    Journal of Bacteriology (2014) Volume 196, Issue 12.
    DOI: https://doi.org/10.1128/JB.01688-14.

    :param source: the source genome to be compared
    :type source: Genome
    :param target: the target genome to be compared
    :type target: Genome
    :return: similarity
    """
    shared_phams = source.intersection(target)

    s_conserved = sum([source[x] for x in shared_phams])
    t_conserved = sum([target[x] for x in shared_phams])
    s_total = len(source)
    t_total = len(target)

    similarity = float(s_conserved + t_conserved) / (s_total + t_total)

    return round(similarity * 100, 2)


def gene_content_similarity(source, target):
    """Calculate the gene content similarity between the given source
    and target genomes.

    Mavrich and Hatfull. Bacteriophage Evolution Differs by Host,
    Lifestyle, and Genome.
    Nature Microbiology (2017) Volume 2.
    DOI: https://doi.org/10.1038/nmicrobiol.2017.112.

    :param source: the source genome to be compared
    :type source: Genome
    :param target: the target genome to be compared
    :type target: Genome
    :return: similarity
    """
    shared_phams = source.intersection(target)

    source_gcs = float(len(shared_phams)) / len(source)
    target_gcs = float(len(shared_phams)) / len(target)

    similarity = average([source_gcs, target_gcs])

    return round(similarity * 100, 2)


def jaccard_similarity_coefficient(source, target):
    """Calculate the jaccard similarity coefficient between the given
    source and target genomes.

    Jaccard. The Distribution of the Flora in the Alpine Zone.
    New Phytologist (1912) Volume 11, Issue 2.
    DOI: https://doi.org/10.1111/j.1469-8137.1912.tb05611.x.

    :param source: the source genome to be compared
    :type source: Genome
    :param target: the target genome to be compared
    :type target: Genome
    :return: similarity
    """
    shared_phams = source.intersection(target)
    all_phams = source.union(target)

    similarity = float(len(shared_phams)) / len(all_phams)

    return round(similarity * 100, 2)
