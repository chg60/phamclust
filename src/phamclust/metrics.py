"""Functions for calculating various similarity metrics between pairs
of genomes.

Alignment-free metrics:

- gene content similarity
- jaccard coefficient
- percentage of conserved proteins
- alignment fraction

Alignment-based metrics:

- average amino acid identity
- proteomic equivalence quotient

All metrics can be optionally returned as distances rather than
similarities.
"""

from parasail import blosum62, nw_trace_diag_16

from phamclust.genome import Genome
from phamclust.statistics import average


def gene_content_similarity(source, target, as_distance=False):
    """Calculate the gene content similarity between the given source
    and target genomes, as described in Mavrich and Hatfull (2017).

    DOI: https://doi.org/10.1038/nmicrobiol.2017.112.

    :param source: the source genome to be compared
    :type source: Genome
    :param target: the target genome to be compared
    :type target: Genome
    :param as_distance: return distance instead of similarity
    :type as_distance: bool
    :return: similarity
    """
    shared_phams = source.intersection(target)

    if not shared_phams:
        similarity = 0.0
    else:
        numerator = 2.0 * len(shared_phams)
        denominator = len(source.phams) + len(target.phams)

        similarity = numerator / denominator

    if as_distance:
        return round(1.0 - similarity, 6)

    return round(similarity, 6)


def jaccard_coefficient(source, target, as_distance=False):
    """Calculate the jaccard similarity coefficient between the given
    source and target genomes, as described in Jaccard (1912).

    DOI: https://doi.org/10.1111/j.1469-8137.1912.tb05611.x.

    :param source: the source genome to be compared
    :type source: Genome
    :param target: the target genome to be compared
    :type target: Genome
    :param as_distance: return distance instead of similarity
    :type as_distance: bool
    :return: similarity
    """
    shared_phams = source.intersection(target)

    if not shared_phams:
        similarity = 0.0
    else:
        similarity = len(shared_phams) / len(source.union(target))

    if as_distance:
        return round(1.0 - similarity, 6)

    return round(similarity, 6)


def percentage_of_conserved_proteins(source, target, as_distance=False):
    """Calculate the percentage of conserved proteins between the
    given source and target genomes, as described in Qin et al. (2014).

    DOI: https://doi.org/10.1128/JB.01688-14.

    :param source: the source genome to be compared
    :type source: Genome
    :param target: the target genome to be compared
    :type target: Genome
    :param as_distance: return distance instead of similarity
    :type as_distance: bool
    :return: similarity
    """
    shared_phams = source.intersection(target)

    if not shared_phams:
        similarity = 0.0
    else:
        s_conserved = sum([len(source[x]) for x in shared_phams])
        t_conserved = sum([len(target[x]) for x in shared_phams])
        numerator = s_conserved + t_conserved

        s_total = len(source)
        t_total = len(target)
        denominator = s_total + t_total

        similarity = numerator / denominator

    if as_distance:
        return round(1.0 - similarity, 6)

    return round(similarity, 6)


def alignment_fraction(source, target, as_distance=False):
    """Calculate the percentage of conserved proteins between the given
    source and target genomes, weighted by gene size.

    :param source: the source genome to be compared
    :type source: Genome
    :param target: the target genome to be compared
    :type target: Genome
    :param as_distance: return distance instead of similarity
    :type as_distance: bool
    :return: similarity
    """
    shared_phams = source.intersection(target)

    if not shared_phams:
        similarity = 0.0
    else:
        s_conserved, s_total = 0, 0
        for pham, translations in source:
            for translation in translations:
                if pham in shared_phams:
                    s_conserved += len(translation)  # tally conserved length
                s_total += len(translation)  # tally total length

        t_conserved, t_total = 0, 0
        for pham, translations in target:
            for translation in translations:
                if pham in shared_phams:
                    t_conserved += len(translation)  # tally conserved length
                t_total += len(translation)  # tally total length

        numerator = s_conserved + t_conserved
        denominator = s_total + t_total

        similarity = numerator / denominator

    if as_distance:
        return round(1.0 - similarity, 6)

    return round(similarity, 6)


def _align(seq_a, seq_b, gap_open=11, gap_extend=1, subst_mat=blosum62):
    """Globally align a pair of sequences.

    :param seq_a: the first sequence to align
    :type seq_a: str
    :param seq_b: the second sequence to align
    :type seq_b: str
    :param gap_open: the gap open penalty
    :type gap_open: int
    :param gap_extend: the gap extend penalty
    :type gap_extend: int
    :param subst_mat: the substitution matrix to use
    :return: traceback
    """
    result = nw_trace_diag_16(seq_a, seq_b, gap_open, gap_extend, subst_mat)
    return result.get_traceback(mch="|", sim="+", neg=" ")


def average_aminoacid_identity(source, target, ppos=False, as_distance=False):
    """Approximate the average amino acid identity (aai) between two
    genomes by computing the length-weighted amino acid identity
    between genes in phams shared by the genomes.

    Use `positive=True` to return the average amino acid similarity
    instead of identity.

    :param source: the source genome to be compared
    :type source: Genome
    :param target: the target genome to be compared
    :type target: Genome
    :param as_distance: return distance instead of similarity
    :type as_distance: bool
    :param ppos: calculate percent positive instead of identical
    :type ppos: bool
    :return: similarity
    """
    shared_phams = source.intersection(target)

    if not shared_phams:
        similarity = 0.0
    else:
        identities, aln_lens = list(), list()

        for pham in shared_phams:
            source_genes = source[pham]
            target_genes = target[pham]

            # Anchor on whichever genome has fewer genes in the pham
            if len(source_genes) > len(target_genes):
                source_genes, target_genes = target_genes, source_genes

            # Find and use the best match for each anchor gene
            for source_gene in source_genes:
                alns = list()
                for target_gene in target_genes:
                    traceback = _align(source_gene, target_gene)
                    length = len(traceback.query)
                    identical_residues = float(traceback.comp.count("|"))
                    # Allow aai to report percent positive if requested
                    if ppos:
                        identical_residues += float(traceback.comp.count("+"))
                    alns.append((identical_residues / length, length))

                best_ident, best_length = sorted(alns, key=lambda x: x[0])[-1]
                identities.append(best_ident)
                aln_lens.append(best_length)

        similarity = average(identities, weights=aln_lens)

    if as_distance:
        return round(1.0 - similarity, 6)

    return round(similarity, 6)


def proteomic_equivalence_quotient(source, target, as_distance=False):
    """Compute the proteomic equivalence quotient between a pair of
    genomes.

    :param source: the source genome to be compared
    :type source: Genome
    :param target: the target genome to be compared
    :type target: Genome
    :param as_distance: return distance instead of similarity
    :type as_distance: bool
    :return: similarity
    """
    similarity = alignment_fraction(source, target)
    similarity *= average_aminoacid_identity(source, target)

    if as_distance:
        return round(1.0 - similarity, 6)

    return round(similarity, 6)


__all__ = ["alignment_fraction", "average_aminoacid_identity",
           "gene_content_similarity", "jaccard_coefficient",
           "percentage_of_conserved_proteins", "proteomic_equivalence_quotient"]
