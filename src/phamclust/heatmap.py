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

    fig.update_layout({'yaxis_nticks': len(matrix),
                       'xaxis_nticks': len(matrix)})

    if not filename:
        fig.show()
    elif filename.suffix == ".html":
        fig.write_html(filename)
    else:
        if len(matrix) > 200:
            fig.write_image(filename, scale=0.125, width=12.5 * len(matrix),
                            height=12.5 * len(matrix), validate=True)
        elif len(matrix) > 100:
            fig.write_image(filename, scale=0.25, width=25 * len(matrix),
                            height=25 * len(matrix), validate=True)
        elif len(matrix) > 50:
            fig.write_image(filename, scale=0.5, width=50 * len(matrix),
                            height=50 * len(matrix), validate=True)
        elif len(matrix) > 25:
            fig.write_image(filename, scale=0.75, width=75 * len(matrix),
                            height=75 * len(matrix), validate=True)
        else:
            fig.write_image(filename)

    return filename
