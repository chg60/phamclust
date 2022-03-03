"""Functions to calculate summary statistics."""


def average(values):
    """Calculate the arithmetic mean of the given values.

    :param values: values to average
    :type values: list[float or int]
    :return: average
    """
    return float(sum(values)) / len(values)


def variance(values, mean=None, sample=True):
    """Calculate the variance of the given values.

    :param values: values to get variance of
    :type values: list[float or int]
    :param mean: the mean of the given values
    :type mean: float
    :param sample: sample or whole population?
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

    :param values: values to get variance of
    :type values: list[float or int]
    :param mean: the mean of the given values
    :type mean: float
    :param sample: sample or whole population?
    :type sample: bool
    :return: std_dev
    """
    return variance(values, mean, sample) ** 0.5


def skewness(values, mean=None, sample=True):
    """Calculate the skewness of the given values.

    :param values: values to get variance of
    :type values: list[float or int]
    :param mean: the mean of the given values
    :type mean: float
    :param sample: sample or whole population?
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
