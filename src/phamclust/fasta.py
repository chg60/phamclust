"""Functions for reading from and writing to FASTA files."""


def read_fasta(filename):
    """Read FASTA records from the indicated file.

    To keep memory utilization low during file parsing, yields records
    one at a time.

    :param filename: path to a file to read FASTA records from
    :type filename: str | pathlib.Path | os.PathLike[str]
    :return: records
    :rtype: typing.Generator[tuple[str, str]]
    """
    with open(filename, "r") as fasta_reader:
        # Read the first line
        line = next(fasta_reader)

        # Iterate until EOF
        while True:
            # Parse the header
            if line[0] != ">":
                raise ValueError(f"records in FASTA files must start with '>'")
            header = line[1:].rstrip()

            # Parse the sequence
            sequence = list()
            for line in fasta_reader:
                if line[0] == ">":
                    break   # ">" tells us we are moving to the next record
                sequence.append(line.rstrip())
            else:           # ">" not encountered means we hit EOF
                line = None

            sequence = "".join(sequence)

            yield header, sequence

            if line is None:
                break


def write_fasta(records, filename, mode="w", wrap=None):
    """Write FASTA records to the indicated file.

    Enable line-wrapping for sequence lines by using `wrap` with a
    positive integer.

    :param records: (header, sequence) tuples to write
    :type records: typing.Iterable[tuple[str, str]]
    :param filename: path to a file to write FASTA records to
    :type filename: str | pathlib.Path | os.PathLike[str]
    :param mode: the mode in which to open `filename`
    :type mode: str
    :param wrap: fixed width at which to wrap sequence lines
    :type wrap: int
    :return: fastq
    """
    with open(filename, mode) as fasta_writer:
        for header, sequence in records:
            fasta_writer.write(f">{header}\n")

            # Don't line-wrap if wrap is None or negative
            if wrap is None or wrap <= 0:
                fasta_writer.write(f"{sequence}\n")
                continue

            # Don't count newlines against the character width
            for i in range(0, len(sequence), wrap):
                subseq = sequence[i: i + wrap]
                fasta_writer.write(subseq + "\n")

    return filename
