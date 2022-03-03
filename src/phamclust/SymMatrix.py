from phamclust.multiprocess import parallelize
from phamclust.statistics import average, skewness, standard_deviation


class SymMatrix:
    """A SymMatrix is an undirected graph for storing pairwise edge
    data.

    This class should NOT be used for directed graphs - implementation
    does not preserve directionality of relationships.
    """
    def __init__(self, names):
        """Constructor for SymMatrix instances.

        :param names: node names for this matrix
        :type names: list[str]
        """
        self._node_names = sorted(names)
        self._matrix = dict()

        # Populate source names for O(1) membership testing
        for name in self._node_names:
            self._matrix[name] = dict()

        self.fill_diagonal(100.0)

    @property
    def node_names(self):
        """Return the node names."""
        # return self._matrix.keys()
        return self._node_names[:]

    @property
    def medoid(self):
        """Return the name of the most central node. The medoid node
        has the highest average similarity to all other nodes."""
        medoid, medoid_sum = None, 0
        for source in self._node_names:
            source_sum = sum(self[source].values())
            if source_sum > medoid_sum:
                medoid = source
                medoid_sum = source_sum

        return medoid

    @property
    def epsilon(self):
        """Return the matrix's skew-adjusted mean edge value."""
        if len(self) == 1:
            return 100.0

        edges = list()
        for _, __, weight in self:
            edges.append(weight)

        mean = average(edges)
        std_dev = standard_deviation(edges, mean)
        skew = skewness(edges, mean)

        return mean + (skew * std_dev)

    def reorder(self, names):
        """Use the incoming names as the row order for this matrix.

        :param names:
        :type names: list[str]
        """
        if len(names) != len(self):
            raise ValueError(f"expected len(names) == {len(self)}")

        for name in names:
            if name not in self:
                raise KeyError(f"node ({name}) not in matrix")

        self._node_names = names

    def iterrows(self):
        """Iterate over the matrix rows."""
        pass

    def extract_submatrix(self, names):
        """Return the sub-graph indicated by `names`.

        :param names: names of the nodes to build a submatrix from
        :type names: list[str]
        :return: submatrix
        """
        submatrix = SymMatrix(sorted(names))
        node_names = submatrix.node_names
        for i, source in enumerate(node_names):
            for target in node_names[i+1:]:
                weight = self.get_cell(source, target)
                submatrix.fill_cell(source, target, weight)

        return submatrix

    def nearest_neighbors(self, source, threshold):
        """Return the source node's nearest neighbors above threshold.

        :param source: name of the source node to search
        :type source: str
        :param threshold: similarity threshold to be neighbors
        :type threshold: float
        :return:
        """
        if source not in self:
            raise KeyError(f"node ({source}) not in matrix")

        neighbors = list()
        for target, weight in self[source].items():
            if weight >= threshold and source != target:
                neighbors.append(target)

        return neighbors

    def fill_diagonal(self, value):
        """Fill the matrix diagonal."""
        for source in self._node_names:
            self.fill_cell(source, source, value)

    def fill_cell(self, source, target, value):
        """Fill a specific cell of the matrix."""
        if source not in self:
            raise KeyError(f"source node ({source}) not in matrix")
        if target not in self:
            raise KeyError(f"target node ({target}) not in matrix")

        if not 0 <= value <= 100:
            raise ValueError(f"value ({value}) not in range(0, 100)")

        self._matrix[source][target] = value
        if source != target:
            self._matrix[target][source] = value

    def get_cell(self, source, target):
        """Return a specific cell of the matrix"""
        if source not in self:
            raise KeyError(f"source node ({source}) not in matrix")
        if target not in self:
            raise KeyError(f"target node ({target}) not in matrix")

        return self[source][target]

    def write_clustalo(self, filepath):
        """Write matrix to a file, as a Clustal Omega identity matrix.

        :param filepath: path to the file to write matrix to
        :type filepath: pathlib.Path
        """
        with open(filepath, "w") as matrix_writer:
            matrix_writer.write(str(self))

    def write_adjacency(self, filepath):
        """Write matrix to a file, as an adjacency matrix.

        :param filepath: path to the file to write matrix to
        :type filepath: pathlib.Path
        """
        with open(filepath, "w") as matrix_writer:
            for source, target, weight in self:
                matrix_writer.write(f"{source}\t{target}\t{weight}\n")

    def __str__(self):
        s = f"{len(self)}\n"
        for source in self._node_names:
            s += f"{source:<24}"
            for target in self._node_names:
                s += " " + str(self.get_cell(source, target))
            s += "\n"

        return s

    def __lt__(self, other):
        if not isinstance(other, SymMatrix):
            raise TypeError(f"other expected 'SymMatrix', not '{type(other)}'")

        return len(self) < len(other)

    def __getitem__(self, item):
        if not isinstance(item, str):
            raise TypeError(f"item expected 'str', not '{type(item)}'")

        return self._matrix[item]

    def __contains__(self, item):
        if not isinstance(item, str):
            raise TypeError(f"item expected 'str', not '{type(item)}'")

        return item in self._matrix

    def __iter__(self):
        """Iterate matrix as an adjacency matrix."""
        for i, source in enumerate(self._node_names):
            for target in self._node_names[i:]:
                yield source, target, self._matrix[source][target]

    def __len__(self):
        return len(self._matrix)


def matrix_from_adjacency(genomes, filepath):
    """Create a symmetric matrix by reading in an adjacency matrix
    file.

    :param genomes: the genomes to build a matrix from
    :type genomes: dict[str, PhamClust.Genome.Genome]
    :param filepath: path to an adjacency matrix file
    :type filepath: pathlib.Path
    :return: matrix
    """
    matrix = SymMatrix(list(genomes.keys()))

    with open(filepath, "r") as adjacency_reader:
        for row in adjacency_reader:
            source, target, value = row.rstrip().split()
            matrix.fill_cell(source, target, float(value))

    return matrix


# Helper function for matrix_de_novo() --> facilitates parallelization
def _calculate_adjacency(genomes, sources, targets, func):
    """Calculate the similarity between source and target genomes using
    `func`.

    :param genomes: all genomes in the dataset
    :type genomes: dict[str, PhamClust.Genome.Genome]
    :param sources: source genome list
    :type sources: list[str]
    :param targets: target genome lists
    :type targets: list[list[str]]
    :param func: similarity function to use for genome comparisons
    :type func: function
    """
    if len(sources) != len(targets):
        raise ValueError(f"expected len(sources) == len(targets)")

    adjacency = list()
    for source_name, target_list in zip(sources, targets):
        source_genome = genomes[source_name]
        for target_name in target_list:
            target_genome = genomes[target_name]
            weight = func(source_genome, target_genome)
            adjacency.append((source_name, target_name, weight))

    return adjacency


def matrix_de_novo(genomes, func, cpus, batch_size=20):
    """Create a symmetric matrix from scratch using `func` to measure
    the similarity between genomes.

    :param genomes: the genomes to build a matrix from
    :type genomes: dict[str, PhamClust.Genome.Genome]
    :param func: similarity function to use for genome comparisons
    :type func: function
    :param cpus: number of CPU cores to use for matrix-building
    :type cpus: int
    :param batch_size: approximate number of genomes per
    :type batch_size: int
    :return: matrix
    """
    matrix = SymMatrix(list(genomes.keys()))
    node_names = matrix.node_names

    num_batches = len(genomes) // batch_size
    while len(genomes) / num_batches > batch_size or num_batches % cpus != 0:
        num_batches += 1

    jobs = list()
    for batch_index in range(num_batches):
        batch_indices = range(len(node_names))[batch_index::num_batches]
        batch_sources = [node_names[x] for x in batch_indices]
        batch_targets = [node_names[x+1:] for x in batch_indices]
        jobs.append((genomes, batch_sources, batch_targets, func))

    for batch_result in parallelize(jobs, cpus, _calculate_adjacency):
        for source, target, value in batch_result:
            matrix.fill_cell(source, target, value)

    return matrix
