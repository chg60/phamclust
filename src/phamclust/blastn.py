import shlex
from subprocess import Popen, PIPE


BLASTN_FIELDS = ['qseqid', 'sseqid', 'nident', 'length', 'qstart', 'qend',
                 'sstart', 'send', 'qlen', 'slen', 'evalue' 'bitscore']


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


if __name__ == "__main__":
    import pathlib
    import sys

    from joblib import Parallel, delayed

    from phamclust.parallel_process import CPUS
    from phamclust.matrix import matrix_from_adjacency, matrix_to_squareform

    indir = pathlib.Path(sys.argv[1]).resolve()

    outdir = indir.parent.joinpath("blastn")
    if not outdir.is_dir():
        outdir.mkdir()

    fnas = [x for x in indir.iterdir() if x.suffix == ".fna"]

    # Initialize parallel runner - disable memory mapping
    runner = Parallel(n_jobs=CPUS, return_as="generator_unordered",
                      max_nbytes=None)

    # Set up and dispatch the parallel jobs
    batch, temp_outs = list(), list()
    for q in fnas:
        outfile = outdir.joinpath(f"{q.stem}.tsv")
        if outfile.is_file():
            temp_outs.append(outfile)
            continue
        batch.append((q, fnas, outdir))
    temp_outs.extend(runner(delayed(blastn_multiple)(*b) for b in batch))

    outfile = outdir.parent.joinpath("blastn_adjacency.tsv")
    if not outfile.is_file():
        blastn_map = dict()
        for temp_out in temp_outs:
            with open(temp_out, "r") as temp_reader:
                for row in temp_reader:
                    source, target, *(data) = row.rstrip().split("\t")
                    for i, datum in enumerate(data):
                        try:
                            data[i] = int(datum)
                        except ValueError:
                            data[i] = float(datum)

                    if source in blastn_map:
                        if target in blastn_map[source]:
                            blastn_map[source][target].append(data)
                        else:
                            blastn_map[source][target] = [data]
                    else:
                        blastn_map[source] = {target: [data]}

        fh = open(outfile, "w")

        nodes = sorted(blastn_map.keys())
        for i, source in enumerate(nodes):
            for target in nodes[i:]:
                source_data = blastn_map[source].get(target, dict())
                target_data = blastn_map[target].get(source, dict())
                if not source_data and not target_data:
                    weight = 0.0
                elif not source_data:
                    numerator = sum([int(x[0]) for x in target_data])
                    denominator = target_data[0][-3]
                    weight = min([numerator/denominator, 1.0])
                elif not target_data:
                    numerator = sum([int(x[0]) for x in source_data])
                    denominator = source_data[0][-3]
                    weight = min([numerator/denominator, 1.0])
                else:
                    numerator = sum([int(x[0]) for x in source_data])
                    numerator += sum([int(x[0]) for x in target_data])
                    denominator = int(source_data[0][-3])
                    denominator += int(target_data[0][-3])
                    weight = min([numerator/denominator, 1.0])

                fh.write(f"{source}\t{target}\t{1.0 - weight:.6f}\n")

        fh.close()

    distmat = matrix_from_adjacency(outfile)
    matrix_to_squareform(distmat,
                         outfile.with_name("blastn_distance_matrix.tsv"),
                         lower_triangle=True)
