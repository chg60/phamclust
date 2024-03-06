"""Functions for working with matrix heatmaps."""

import plotly.express as px


def draw_heatmap(matrix, midpoint=None, colorblind=False, filename=None):
    """Create a heatmap representation of a matrix.

    :param matrix: the matrix to represent as a heatmap
    :type matrix: phamclust.matrix.SymMatrix
    :param midpoint:
        Use a 3-color (red/white/green) scale with white at `midpoint`.
        If `midpoint` is None, heatmap will be grayscale.
    :type midpoint: float
    :param colorblind: use colorblind-friendlier colors in heatmap
    :type colorblind: bool
    :param filename: name of the file to save the heatmap to
    :type filename: pathlib.Path
    :return: filename
    """
    labels = matrix.nodes

    if not midpoint:
        scale = [(0, "white"), (1, "black")]
    elif colorblind:
        scale = [(0, "magenta"), (midpoint, "white"), (1, "turquoise")]
    else:
        scale = [(0, "red"), (midpoint, "yellow"), (1, "green")]

    # noinspection PyTypeChecker
    fig = px.imshow(matrix.to_ndarray(), x=labels, y=labels,
                    color_continuous_scale=scale, range_color=[0.0, 1.0])

    # Values have been empirically derived using matrices from size 2-750
    width, height = 16 * len(matrix), 9 * len(matrix)
    scale = 960 / height

    while scale > 5:
        width, height = width * 2, height * 2
        scale = 960 / height

    font_size = int(max([8.5 - 2 * scale, 5.67 - (2 / 3 * scale)]))

    fig.update_layout({'yaxis_nticks': len(matrix),
                       'xaxis_nticks': len(matrix),
                       'font_size': max([1, font_size])})

    if not filename:
        fig.show()
    elif filename.suffix == ".html":
        # Set canvas size based on screen size for HTML output
        fig.write_html(filename, default_height='100%', default_width='100%')
    else:
        # Use calculated scale, width, height for all others
        fig.write_image(filename, scale=scale, width=width, height=height,
                        validate=True)

    return filename
