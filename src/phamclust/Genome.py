"""A dictionary-based class for efficiently storing the gene content
encoded by genomes, and intuitively comparing genomes to each other."""


class Genome:
    """A Genome is a collection of gene phamilies that may or may not
    be shared by other Genomes. Paralogs (multiple copies of the same
    phamily) are tolerated.

    Class methods allow this object to behave like something between a
    set and a dictionary."""
    def __init__(self, name, phams=None):
        """Constructor for Genome instances.

        :param name: name for this genome
        :type name: str
        :param phams: phamilies found in this genome
        :type phams: list[str]
        """
        self.name = name
        self.phams = dict()

        if phams:
            if not isinstance(phams, list):
                raise TypeError(f"phams expected 'list', not '{type(phams)}'")

            for pham in phams:
                self.add(pham)

    def add(self, pham):
        """Add a pham to this genome.

        :param pham: name of the pham to add to this genome
        :type pham: str
        """
        if not isinstance(pham, str):
            raise TypeError(f"pham expected 'str', not '{type(pham)}'")

        if pham not in self:
            self.phams[pham] = 1
        else:
            self.phams[pham] += 1

    def pop(self, pham):
        """Pop a pham from this genome.

        :param pham: name of the pham to pop from this genome
        :type pham: str
        """
        if not isinstance(pham, str):
            raise TypeError(f"pham expected 'str', not '{type(pham)}'")

        self.phams.pop(pham)

    def intersection(self, other):
        """Get the set of phams present both in this Genome and
        another Genome.

        :param other: another genome to be checked for shared phams
        :type other: Genome
        """
        if not isinstance(other, Genome):
            raise TypeError(f"other expected 'Genome', not '{type(other)}'")

        return set(self.phams.keys()).intersection(other.phams.keys())

    def union(self, other):
        """Get the set of phams present in either this Genome or
        another Genome.

        :param other: another genome to be checked for shared phams
        :type other: Genome
        """
        if not isinstance(other, Genome):
            raise TypeError(f"other expected 'Genome', not '{type(other)}'")

        return set(self.phams.keys()).union(other.phams.keys())

    def difference(self, other):
        """Get the set of phams present in this Genome but not another
        Genome.

        :param other: another genome to be checked for shared phams
        :type other: Genome
        """
        if not isinstance(other, Genome):
            raise TypeError(f"other expected 'Genome', not '{type(other)}'")

        return set(self.phams.keys()).difference(other.phams.keys())

    def symmetric_difference(self, other):
        """Get the set of phams present in either this Genome or
        another, but not both.

        :param other: another genome to be checked for shared phams
        :type other: Genome
        """
        if not isinstance(other, Genome):
            raise TypeError(f"other expected 'Genome', not '{type(other)}'")

        return set(self.phams.keys()).symmetric_difference(other.phams.keys())

    def __lt__(self, other):
        if not isinstance(other, Genome):
            raise TypeError(f"other expected 'Genome', not '{type(other)}'")

        return len(self) < len(other)

    def __getitem__(self, item):
        if not isinstance(item, str):
            raise TypeError(f"item expected 'str', not '{type(item)}'")

        return self.phams[item]

    def __contains__(self, item):
        if not isinstance(item, str):
            raise TypeError(f"item expected 'str', not '{type(item)}'")

        return item in self.phams

    def __iter__(self):
        for key, value in self.phams.items():
            yield key, value

    def __len__(self):
        return sum([value for key, value in self])
