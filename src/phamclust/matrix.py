"""Symmetric matrix class to store pairwise relational data.

Includes some helper functions to allow loading SymMatrix objects from
files or saving them to files. The following representations are
supported: adjacency, squareform, or lower triangle matrices.

Adjacency matrices should be 3-column tables: source, target, weight.
The source and target fields are interpreted as strings, while the
weight is cast to a float. This is the least dense but most portable
way to represent a symmetric matrix.

Squareform matrices are denser than adjacency files, but because
symmetric matrices store redundant information (matrix[source,target]
is the same as matrix[target,source]), this is not the most efficient
representation.

Lower triangle files are the densest matrix representation supported
here. They include the matrix diagonal to allow inference of whether
the matrix is a distance or similarity matrix.
"""

# from itertools import combinations_with_replacement
import logging

import joblib
from numpy import array

# from phamclust.parallel_process import parallelize
from phamclust.statistics import average, standard_deviation, skewness


class SymMatrix:
    """An undirected graph for storing pairwise edge weights.

    Can be used for identity/similarity matrices or distance matrices.

    Should NOT be used for directed graphs, as edge directionality is
    not retained.
    """
    def __init__(self, nodes, is_distance=False):
        """Constructor method.

        :param nodes: node names to be used for this matrix
        :type nodes: list[str]
        """
        self._nodes = nodes
        self._matrix = {node: dict() for node in nodes}
        self._is_distance = is_distance
        self._locked = False

    @property
    def is_distance(self):
        """Return whether this is a distance matrix."""
        return self._is_distance

    @property
    def nodes(self):
        """Return a deep copy of the node names."""
        return self._nodes[:]

    @property
    def diameter(self):
        """Return the cluster diameter (the longest edge weight).

        Because diameter is inherently a measure of distance, a
        ValueError will be raised if invoked for a similarity matrix.
        """
        if not self.is_distance:
            raise ValueError(f"cannot compute diameter for similarity matrix")

        diameter = 0.0

        for source, target, weight in self:
            if source != target and weight > diameter:
                diameter = weight

        return round(diameter, 6)

    @property
    def medoid(self):
        """Return the name of the most central node, and its average
        similarity/distance to all other nodes."""
        temp_dict = dict()
        for source, target, weight in self:
            if source in temp_dict:
                temp_dict[source].append(weight)
            else:
                temp_dict[source] = [weight]

            if target in temp_dict:
                temp_dict[target].append(weight)
            else:
                temp_dict[target] = [weight]

        for node, weights in temp_dict.items():
            temp_dict[node] = average(weights)

        if self.is_distance:
            medoid = sorted(temp_dict.items(), key=lambda x: x[1])[0]
        else:
            medoid = sorted(temp_dict.items(), key=lambda x: x[1])[-1]

        return medoid

    @property
    def anti_medoid(self):
        """Return the name of the least central node, and its average
        similarity/distance to all other nodes."""
        temp_dict = dict()
        for source, target, weight in self:
            if source in temp_dict:
                temp_dict[source].append(weight)
            else:
                temp_dict[source] = [weight]

            if target in temp_dict:
                temp_dict[target].append(weight)
            else:
                temp_dict[target] = [weight]

        for node, weights in temp_dict.items():
            temp_dict[node] = average(weights)

        if self.is_distance:
            anti = sorted(temp_dict.items(), key=lambda x: x[1])[-1]
        else:
            anti = sorted(temp_dict.items(), key=lambda x: x[1])[0]

        return anti

    @property
    def statistics(self):
        """Return the mean, standard deviation, skewness, and median of
        the unique edges of this matrix."""
        if len(self) == 1:
            node = self._nodes[0]
            return self.get_weight(node, node), 0.0, 0.0

        # Collect unique edges, skipping the diagonal (self edges)
        edges = list()
        for source, target, weight in self:
            if source != target:
                edges.append(weight)

        mean = average(edges)
        std_dev = standard_deviation(edges, mean)
        if std_dev == 0.0:
            skew = 0.0
        else:
            skew = skewness(edges, mean)

        return mean, std_dev, skew

    def extract_submatrix(self, nodes):
        """Return the sub-graph indicated by `nodes`.

        :param nodes: nodes from which to build a new submatrix
        :type nodes: list[str]
        """
        submatrix = SymMatrix(nodes, self.is_distance)
        for i, source in enumerate(nodes):
            for target in nodes[i:]:
                weight = self.get_weight(source, target)
                submatrix.set_weight(source, target, weight)

        return submatrix

    def append_node(self, source, data):
        """Append a new node to the end of this matrix.

        Caller must have pre-calculated the similarity/distance of the
        new node to all existing nodes.

        Integrity checks:

        * new name is NOT already present in matrix node names
        * ALL current node names are present in incoming data
        * `data[name]` should match the matrix diagonal (0 or 1)

        :param source: the new node's name
        :type source: str
        :param data: the new node's similarity to existing nodes
        :type data: dict[str, float]
        """
        if source in self:
            raise KeyError(f"node '{source}' is already in this matrix")

        if source not in data:
            raise KeyError(f"incoming data lacks an self-edge for '{source}'")

        if self.is_distance and data[source] != 0.0:
            raise ValueError(f"nonsense value {data[source]} for self-edge on "
                             f"distance matrix")
        elif not self.is_distance and data[source] != 1.0:
            raise ValueError(f"nonsense value {data[source]} for self-edge on "
                             f"similarity matrix")

        missed_nodes = set(self._matrix.keys()) - data.keys()
        if len(missed_nodes) > 0:
            raise KeyError(f"missing edge(s) for {source} vs: {missed_nodes}")

        extra_nodes = set(data.keys()) - {source} - set(self._matrix.keys())
        if len(extra_nodes) > 0:
            raise KeyError(f"specified edges for nodes not found in matrix: "
                           f"{extra_nodes}")

        # If we got here, the data are probably safe to use
        self._nodes.append(source)
        self._matrix[source] = dict()

        for target in self._nodes:
            self.set_weight(source, target, data[target])

    def get_weight(self, source, target):
        """Return the weight stored between a pair of nodes.

        :param source: the source node name
        :type source: str
        :param target: the target node name
        :type target: str
        :return: weight stored between source and target
        :rtype: float
        """
        if source not in self:
            raise KeyError(f"node '{source}' not in matrix")
        elif target not in self:
            raise KeyError(f"node '{target}' not in matrix")

        # use lexicographically smaller node as source, larger as target
        if source > target:
            source, target = target, source

        return self._matrix[source].get(target, None)

    def invert(self):
        """Flip the matrix from a distance to identity matrix (or
        vice-versa), in-place.

        Iterates over all edges subtracting their weights from 1.
        """
        for source, target, weight in self:
            self.set_weight(source, target, 1.0 - weight)

        self._is_distance = not self.is_distance

        return self

    def reorder(self, nodes):
        """Change the matrix traversal order, using the incoming nodes."""
        # Sanity checks: every node in matrix should also be in nodes
        if len(nodes) != len(self):
            raise ValueError(f"need {len(self)} nodes but got {len(nodes)}")

        for node in nodes:
            if node not in self:
                raise KeyError(f"node '{node}' not in matrix")

        self._nodes = nodes

    def nearest_neighbors(self, source, threshold):
        """Return the source node's nearest neighbors closer than
        threshold.

        Threshold is interpreted as similarity if this is not a
        distance matrix. Otherwise, threshold is interpreted as a
        distance.

        :param source: name of the source node to search
        :type source: str
        :param threshold: similarity threshold to be neighbors
        :type threshold: float
        :return: neighbors
        """
        is_distance = self.is_distance

        neighbors = list()
        for target in self._nodes:
            if source == target:
                continue

            weight = self.get_weight(source, target)

            if is_distance and weight <= threshold:
                neighbors.append(target)
            elif not is_distance and weight >= threshold:
                neighbors.append(target)

        # Sort in descending order of similarity / ascending order of distance
        return sorted(neighbors, reverse=not is_distance,
                      key=lambda x: self.get_weight(source, x))

    def set_weight(self, source, target, weight):
        """Set the weight between a pair of nodes.

        :param source: the source node name
        :type source: str
        :param target: the target node name
        :type target: str
        :param weight: the edge weight between source and target
        :type weight: float
        :rtype: None
        """
        if self._locked:
            raise AttributeError("matrix is marked as read-only")

        if source not in self:
            raise KeyError(f"node '{source}' not in matrix")
        elif target not in self:
            raise KeyError(f"node '{target}' not in matrix")

        if not 0 <= weight <= 1:
            raise ValueError(f"weight {weight} not in [0.0, 1.0]")

        # use lexicographically smaller node as source, larger as target
        if source > target:
            source, target = target, source

        self._matrix[source][target] = round(weight, 6)

    def lock(self):
        """Mark this matrix as read-only."""
        self._locked = True

    def unlock(self):
        """Mark this matrix as read-write."""
        self._locked = False

    def is_locked(self):
        """Check whether this matrix is read-only."""
        return self._locked

    def to_ndarray(self):
        """Cast this matrix to a numpy array."""
        return array([x for _, x in self.iterrows()])

    def __contains__(self, item):
        """Check whether a node is found in the matrix."""
        if not isinstance(item, str):
            raise TypeError(f"type(item) should be 'str', not '{type(item)}'")

        return item in self._matrix

    def __getitem__(self, item):
        """Accessor for matrix rows.

        :param item: a node to retrieve matrix data for
        :type item: str
        :rtype: dict[str, float]
        """
        if item not in self:
            raise KeyError(f"node '{item}' not in matrix")

        return self._matrix[item]

    def __iter__(self):
        """Iterate over the matrix as an adjacency matrix.

        NOTE: only returns the upper/right half.

        :rtype: typing.Generator[tuple[str, str, float]]
        """
        # Loop over all source nodes
        for idx, source in enumerate(self._nodes):
            # Loop over all target nodes at or after the source
            for target in self._nodes[idx:]:
                yield source, target, self.get_weight(source, target)

    def iterrows(self):
        """Iterate over matrix edges.

        :rtype: typing.Generator[tuple[str, list[float]]]"""
        for source in self._nodes:
            row = [self.get_weight(source, target) for target in self._nodes]
            yield source, row

    def __len__(self):
        """Return the number of nodes in the matrix."""
        return len(self._nodes)

    def __lt__(self, other):
        if not isinstance(other, SymMatrix):
            raise TypeError(f"cannot compare SymMatrix to {type(other)}")

        return len(self) < len(other)

    def __str__(self):
        s = f"{len(self)}\n"
        for source, row in self.iterrows():
            s += f"{source:<24}\t"
            s += "\t".join([f"{x:.6f}" for x in row]) + "\n"

        return s


# SymMatrix de novo
def _outside_in_index_iterator(num_indices):
    """Generator that iterates over list indices from the front and back
    toward the middle.

    :param num_indices: the size of the list to be iterated
    :type num_indices: int
    >>> list(_outside_in_index_iterator(5))
    [0, 4, 1, 3, 2]
    """
    start, mid, end = 0, num_indices // 2, num_indices - 1
    for x, y in zip(range(start, mid), range(end, mid - 1, -1)):
        yield x
        yield y
    if num_indices % 2 == 1:
        yield mid


def calculate_adjacency(source, target, func, distance=True):
    """Calculate the distance between two genomes using the provided
    function."""
    return source.name, target.name, func(source, target, as_distance=distance)


def matrix_de_novo(genomes, func, cpus, as_distance=True):
    """Create a symmetric matrix from scratch using `func` to measure
    the similarity between genomes.

    Developer note: this list[Genome]-based implementation has similar
    memory footprint to the previous dict[str, Genome]-based approach,
    but is ~25-30% faster.

    :param genomes: the genomes to build a matrix from
    :type genomes: list[phamclust.genome.Genome]
    :param func: similarity function to use for genome comparisons
    :type func: function
    :param cpus: number of CPU cores to use for matrix-building
    :type cpus: int
    :param as_distance: populate matrix with distance, not similarity
    :type as_distance: bool
    :return: matrix
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

    # <<< Initialize distance matrix, populate the diagonal for "free" >>>
    matrix = SymMatrix(nodes=[g.name for g in genomes], is_distance=as_distance)
    for genome in genomes:
        matrix.set_weight(genome.name, genome.name, 1.0 - as_distance)

    # Initialize parallel runner - disable memory mapping
    runner = joblib.Parallel(n_jobs=cpus, return_as="generator",
                             max_nbytes=None)

    # TODO: a future release of joblib - generator_unordered will be faster!!
    # runner = joblib.Parallel(n_jobs=cpus, return_as="generator_unordered",
    #                          max_nbytes=None)

    # Determine how many batches to chunk the data into - target
    # 10,000 calculations per cpu per batch
    stride = int(len(genomes) // (n_pairs / (10000 * cpus)))
    iter_order = [i for i in _outside_in_index_iterator(len(genomes))]
    batches = [iter_order[i:i+stride] for i in range(0, len(genomes), stride)]
    for batch_indices in batches:
        sources = [genomes[i] for i in batch_indices]
        targets = [genomes[i+1:] for i in batch_indices]

        batch = list()
        for source, targets in zip(sources, targets):
            for target in targets:
                batch.append((source, target, func, as_distance))

        results = runner(joblib.delayed(calculate_adjacency)(*b) for b in batch)

        for _source, _target, _weight in results:
            matrix.set_weight(_source, _target, _weight)

        logging.debug(f"finished {', '.join([x.name for x in sources])}")

    del runner

    return matrix


# Adjacency matrix I/O
def matrix_from_adjacency(filepath):
    """Construct a SymMatrix from an adjacency matrix in a file.

    This is somewhat inefficient because the SymMatrix constructor
    needs to know what the node names are up front. As a result, this
    function streams the file 1 + (1/N) times - stopping once it finds
    a source node that's already been seen as a target node. The second
    pass reads the whole file and populates the pairwise edges.

    This approach gets around the need to allocate a potentially large
    chunk of memory to store and iterate the full adjacency matrix
    before building the SymMatrix.

    :param filepath: the file to read the adjacency matrix from
    :type filepath: pathlib.Path
    :return: matrix
    """
    # Learn node names and whether this is a distance matrix
    names, diagonal = dict(), set()
    for source, target, weight in read_adjacency(filepath):
        if source == target:
            diagonal.add(weight)

        if target in names:
            break

        names[target] = None

    names = list(names.keys())

    # Analyze diagonal: all values should be the same
    if len(diagonal) > 1:
        raise ValueError(f"values on matrix diagonal should be identical")

    # Analyze diagonal: should be either 0.0 (distance) or 1.0 (similarity)
    if sum(diagonal) == 0.0:
        is_distance = True
    elif sum(diagonal) == 1.0:
        is_distance = False
    else:
        raise ValueError(f"values on matrix diagonal can only be 0.0 or 1.0")

    # Create the matrix
    matrix = SymMatrix(names, is_distance=is_distance)
    for source, target, weight in read_adjacency(filepath):
        matrix.set_weight(source, target, weight)

    return matrix


def read_adjacency(filepath):
    """Read edges from an adjacency matrix in a file.

    :param filepath: the file to read the adjacency matrix from
    :type filepath: pathlib.Path
    :rtype: typing.Generator[tuple[str, str, float]]
    """
    with open(filepath, "r") as adjacency_reader:
        for edge in adjacency_reader:
            edge = edge.rstrip().split("\t")
            yield edge[0], edge[1], float(edge[2])


def matrix_to_adjacency(matrix, filepath):
    """Write matrix to a file, as an adjacency matrix.

    :param matrix: the matrix to write to a file
    :type matrix: SymMatrix
    :param filepath: path to the file to write matrix to
    :type filepath: pathlib.Path
    :return: filepath
    """
    with open(filepath, "w") as adjacency_writer:
        for source, target, weight in matrix:
            adjacency_writer.write(f"{source}\t{target}\t{weight:.6f}\n")

    return filepath


# Squareform and lower-triangle matrix I/O
def matrix_from_squareform(filepath):
    """Construct a SymMatrix from a Phylip squareform matrix.

    This is somewhat inefficient because the SymMatrix constructor
    needs to know what the node names are up front. As a result, this
    function streams the squareform matrix file twice. The first pass
    only records the unique node names and is used to initialize the
    SymMatrix. The second pass populates the pairwise edges.

    This approach gets around the need to allocate a potentially large
    chunk of memory to store the full squareform matrix before building
    the SymMatrix.

    :param filepath: the file to read the squareform matrix from
    :type filepath: pathlib.Path
    :return: matrix
    """
    names, diagonal = list(), set()
    for i, (target, row) in enumerate(read_squareform(filepath)):
        names.append(target)
        diagonal.add(row[i])

    # Analyze diagonal: all values should be the same
    if len(diagonal) > 1:
        raise ValueError(f"values on matrix diagonal should be identical")

    # Analyze diagonal: should be either 0.0 (distance) or 1.0 (similarity)
    if sum(diagonal) == 0.0:
        is_distance = True
    elif sum(diagonal) == 1.0:
        is_distance = False
    else:
        raise ValueError(f"values on matrix diagonal can only be 0.0 or 1.0")

    # Create the matrix
    matrix = SymMatrix(names, is_distance=is_distance)
    for i, (target, row) in enumerate(read_squareform(filepath)):
        for (source, weight) in zip(names[:i + 1], row[:i + 1]):
            matrix.set_weight(source, target, weight)

    return matrix


def read_squareform(filepath):
    """Read rows from a squareform matrix in a file.

    :param filepath: the file to read the squareform matrix from
    :type filepath: pathlib.Path
    :rtype: typing.Generator[tuple[str, list[float]]]
    """
    with open(filepath, "r") as squareform_reader:
        next(squareform_reader)     # skip first line, only has matrix size
        for row in squareform_reader:
            row = row.rstrip().split()
            name = row[0]
            row = [float(x) for x in row[1:]]
            yield name, row


def matrix_to_squareform(matrix, filepath, lower_triangle=False):
    """Write matrix to a file, as a squareform matrix.

    :param matrix: the matrix to write to a file
    :type matrix: SymMatrix
    :param filepath: the file to write matrix to
    :type filepath: pathlib.Path
    :param lower_triangle: only write the lower triangle
    :type lower_triangle: bool
    :return: filepath
    """
    with open(filepath, "w") as squareform_writer:
        squareform_writer.write(f"{len(matrix)}\n")

        for i, (source, row) in enumerate(matrix.iterrows()):
            if not lower_triangle:
                row = "\t".join([f"{x:.6f}" for x in row])
            else:
                row = "\t".join([f"{x:.6f}" for x in row[:i + 1]])

            squareform_writer.write(f"{source:<24}\t{row}\n")

    return filepath
