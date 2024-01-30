"""A dict-based class for efficiently storing the gene content encoded
by genomes, and intuitively comparing genomes to each other."""

from phamclust.fasta import read_fasta


class GenomeLoadError(Exception):
    """Error raised when a Genome could not be loaded."""
    pass


class Genome:
    """A Genome is a collection of gene phams which may or may not
    be shared by other Genomes.

    A genome may encode more than one protein belonging to the same
    pham.

    This implementation behaves like a cross between a dictionary and
    a set.
    """
    def __init__(self, name):
        """Constructor method.

        :param name: a name to identify this genome
        :type name: str
        """
        self.name = name
        self.phams = dict()

    def add(self, pham, translation="M"):
        """Add a pham to this genome.

        :param pham: the name of a pham to add to this genome
        :type pham: str
        :param translation: sequence of the gene being added
        :type translation: str
        """
        if not isinstance(pham, str):
            raise TypeError(f"type(pham) should be 'str', not '{type(pham)}'")
        elif not isinstance(translation, str):
            raise TypeError(f"type(translation) should be 'str', not "
                            f"'{type(translation)}'")

        if pham not in self:
            self.phams[pham] = [translation]
        else:
            self.phams[pham].append(translation)

    def difference(self, other):
        """Identify the set of phams found in only `self`, not `other`.

        :param other: a genome to compare to this one
        :type other: Genome
        """
        return self - other

    def intersection(self, other):
        """Identify the set of phams found in both `self` and `other`.

        :param other: a genome to compare to this one
        :type other: Genome
        """
        return self & other

    def load(self, fasta):
        """Load a genome from a FASTA file.

        :param fasta: the path to a single-genome FASTA file
        :type fasta: str | pathlib.Path | os.PathLike[str]
        """
        for header, translation in read_fasta(fasta):
            pham = None
            for field in header.split("|"):
                key, value = field.split("=")
                if key == "pham":
                    pham = value
                    break

            if pham is None:
                raise GenomeLoadError(f"unable to get pham from FASTA header")

            self.add(pham, translation)

    def pop(self, pham):
        """Remove and return a pham from this genome.

        :param pham: name of the pham to remove
        :type pham: str
        """
        if pham not in self:
            return None

        return self.phams.pop(pham)

    def save(self, filepath):
        """Save this genome to a file.

        :param filepath: path to which this genome should be written
        :type filepath: str | pathlib.Path | os.PathLike[str]
        :return:
        """
        with open(filepath, "w") as genome_saver:
            genome_saver.write(str(self))

        return filepath

    def symmetric_difference(self, other):
        """Identify the set of phams found in only `self` or `other`.

        :param other: a genome to compare to this one
        :type other: Genome
        """
        return self ^ other

    def union(self, other):
        """Identify the set of phams found in either `self` or `other`.

        :param other: a genome to compare to this one
        :type other: Genome
        """
        return self | other

    def __and__(self, other):
        """Overload '&' operator to compute the intersection between
        genomes.

        :param other: a genome to compare to this one
        :type other: Genome
        :rtype: set[str]
        """
        if not isinstance(other, Genome):
            raise TypeError(f"cannot compare Genome to '{type(other)}'")

        return set(self.phams.keys()).intersection(other.phams.keys())

    def __contains__(self, item):
        """Test this genome for membership of a particular pham.

        :param item: a pham to challenge the genome with
        :type item: str
        """
        if not isinstance(item, str):
            raise TypeError(f"type(item) should be 'str', not '{type(item)}'")

        return item in self.phams

    def __getitem__(self, item):
        """Return this genome's data associated with a particular pham.

        :param item: a pham whose data should be retrieved
        :type item: str
        :rtype: list[str]
        """
        if item not in self:
            raise KeyError(f"node '{item}' not in matrix")

        return self.phams[item]

    def __iter__(self):
        """Iterate over genome translations.

        :rtype: typing.Generator[tuple[str,list[str]]]
        """
        for key, value in self.phams.items():
            yield key, value

    def __len__(self):
        return sum([len(value) for _, value in self])

    def __lt__(self, other):
        if not isinstance(other, Genome):
            raise TypeError(f"cannot compare Genome to '{type(other)}'")

        return len(self) < len(other)

    def __or__(self, other):
        """Overload '|' operator to compute the union between genomes.

        :param other: a genome to compare to this one
        :type other: Genome
        :rtype: set[str]
        """
        if not isinstance(other, Genome):
            raise TypeError(f"cannot compare Genome to '{type(other)}'")

        return set(self.phams.keys()).union(other.phams.keys())

    def __repr__(self):
        return str(self)

    def __str__(self):
        """Return a FASTA representation of this genome."""
        s = ""
        for pham, translations in self:
            for i, translation in enumerate(translations):
                s += f">name={self.name}|pham={pham}|n={i+1}\n{translation}\n"

        return s

    def __sub__(self, other):
        """Overload '-' operator to compute the difference between
        genomes.

        :param other: a genome to compare to this one
        :type other: Genome
        :rtype: set[str]
        """
        if not isinstance(other, Genome):
            raise TypeError(f"cannot compare Genome to '{type(other)}'")

        return set(self.phams.keys()).difference(other.phams.keys())

    def __xor__(self, other):
        """Overload '^' operator to compute the symmetric difference
        between genomes.

        :param other: a genome to compare to this one
        :type other: Genome
        :rtype: set[str]
        """
        if not isinstance(other, Genome):
            raise TypeError(f"cannot compare Genome to '{type(other)}'")

        return set(self.phams.keys()).symmetric_difference(other.phams.keys())
