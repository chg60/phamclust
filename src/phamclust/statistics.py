"""Functions to calculate summary statistics."""


def average(values, weights=None):
    """Calculate the weighted arithmetic mean of the given values.

    If weights is not given, this is the same as uniformly weighting
    each value.

    :param values: values to average
    :type values: list[float or int]
    :param weights: weights to apply to the given values
    :type weights: list[float or int]
    :return average
    """
    if not weights:
        weights = [1] * len(values)

    if len(values) != len(weights):
        raise ValueError(f"got {len(values)} values and {len(weights)} weights")

    return float(sum([x * y for x, y in zip(values, weights)])) / sum(weights)


def variance(values, mean=None, sample=True):
    """Calculate the variance of the given values.

    Use `sample`=False to calculate the population variance instead
    of sample variance.

    :param values: values for which to calculate the variance
    :type values: list[float or int]
    :param mean: the mean of the given values
    :type mean: float
    :param sample: calculate the sample variance
    :type sample: bool
    :return: variance
    """
    if len(values) == 1:
        return 0.0

    if not mean:
        mean = average(values)

    numerator = sum([(x - mean) ** 2 for x in values])

    if sample:
        denominator = len(values) - 1
    else:
        denominator = len(values)

    return numerator / denominator


def standard_deviation(values, mean=None, sample=True):
    """Calculate the standard deviation of the given values.

    Use `sample`=False to calculate the population standard deviation
    instead of sample standard deviation.

    :param values: values for which to calculate the standard deviation
    :type values: list[float or int]
    :param mean: the mean of the given values
    :type mean: float
    :param sample: calculate the sample standard deviation
    :type sample: bool
    :return: standard_deviation
    """
    return variance(values, mean, sample) ** 0.5


def skewness(values, mean=None, sample=True):
    """Calculate the skewness of the given values.

    User `sample`=False to calculate the population skewness instead
    of sample skewness.

    :param values: values for which to calculate the skewness
    :type values: list[float or int]
    :param mean: the mean of the given values
    :type mean: float
    :param sample: calculate the sample skewness
    :type sample: bool
    :return: skewness
    """
    if len(values) == 1:
        return 0.0

    if not mean:
        mean = average(values)

    numerator = sum([(x - mean) ** 3 for x in values])
    denominator = standard_deviation(values, mean, sample) ** 3

    if sample:
        denominator *= (len(values) - 1)
    else:
        denominator *= len(values)

    return numerator / denominator
