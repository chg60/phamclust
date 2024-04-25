"""Functions for working with matrix heatmaps."""

import plotly.express as px

CSS_COLORS = {'turquoise', 'dimgrey', 'darkolivegreen', 'darkgoldenrod',
              'tomato', 'dodgerblue', 'mediumblue', 'mediumaquamarine',
              'lightsteelblue', 'cyan', 'azure', 'seashell', 'burlywood',
              'whitesmoke', 'teal', 'antiquewhite', 'deepskyblue',
              'mediumvioletred', 'springgreen', 'palegoldenrod', 'slategrey',
              'lightskyblue', 'wheat', 'khaki', 'silver',
              'lightgoldenrodyellow', 'navajowhite', 'yellowgreen',
              'chocolate', 'greenyellow', 'yellow', 'snow', 'orange',
              'mediumturquoise', 'saddlebrown', 'orchid', 'pink', 'magenta',
              'darkcyan', 'olivedrab', 'white', 'lime', 'aqua', 'red',
              'darkmagenta', 'blueviolet', 'beige', 'crimson', 'darkslategray',
              'goldenrod', 'lemonchiffon', 'sandybrown', 'grey', 'darkviolet',
              'forestgreen', 'cornsilk', 'palegreen', 'midnightblue',
              'rosybrown', 'darksalmon', 'indigo', 'tan', 'slateblue', 'blue',
              'mediumseagreen', 'lightpink', 'ghostwhite', 'papayawhip',
              'coral', 'linen', 'plum', 'deeppink', 'lightgrey',
              'lightslategray', 'oldlace', 'floralwhite', 'mintcream',
              'orangered', 'aliceblue', 'olive', 'cornflowerblue', 'moccasin',
              'ivory', 'darkslategrey', 'cadetblue', 'thistle', 'bisque',
              'mediumslateblue', 'navy', 'mediumspringgreen', 'green',
              'lightseagreen', 'limegreen', 'darkgray', 'darkorange',
              'aquamarine', 'hotpink', 'slategray', 'sienna', 'brown',
              'mediumpurple', 'chartreuse', 'lightgray', 'paleturquoise',
              'mistyrose', 'lightslategrey', 'palevioletred', 'powderblue',
              'lightgreen', 'lawngreen', 'lightcyan', 'lavenderblush',
              'mediumorchid', 'lightyellow', 'darkseagreen', 'skyblue',
              'seagreen', 'lightcoral', 'darkred', 'lavender', 'black',
              'peachpuff', 'honeydew', 'lightsalmon', 'maroon', 'dimgray',
              'darkgrey', 'gainsboro', 'darkslateblue', 'darkturquoise',
              'steelblue', 'darkorchid', 'darkkhaki', 'blanchedalmond',
              'peru', 'indianred', 'salmon', 'darkgreen', 'violet', 'purple',
              'fuchsia', 'lightblue', 'gold', 'royalblue', 'firebrick', 'gray',
              'darkblue'}


def draw_heatmap(matrix, colors=None, midpoint=0.5, filename=None):
    """Create a heatmap representation of a matrix.

    :param matrix: the matrix to represent as a heatmap
    :type matrix: phamclust.matrix.SymMatrix
    :param colors: the colors to use for the heatmap
    :type colors: list[str]
    :param midpoint: the midpoint of the color scale
    :type midpoint: float
    :param filename: name of the file to save the heatmap to
    :type filename: pathlib.Path
    :return: filename
    """
    labels = matrix.nodes
    if colors is None:
        colors = ["red", "yellow", "green"]

    if len(colors) == 2:
        color_scale = [(0, colors[0]), (1, colors[1])]
    elif len(colors) == 3:
        color_scale = [(0, colors[0]), (midpoint, colors[1]), (1, colors[2])]
    else:
        raise ValueError("`colors` parameter must be a list of 2 or 3 color "
                         "names (e.g., ['white', 'green'] or ['red', "
                         "'yellow', 'green'])")

    # noinspection PyTypeChecker
    fig = px.imshow(matrix.to_ndarray(), x=labels, y=labels,
                    color_continuous_scale=color_scale, range_color=[0.0, 1.0])

    # Values have been empirically derived using matrices from size 2-750
    width, height = 16 * len(matrix), 9 * len(matrix)
    px_scale = 960 / height

    while px_scale > 5:
        width, height = width * 2, height * 2
        px_scale = 960 / height

    font_size = int(max([8.5 - 2 * px_scale, 5.67 - (2 / 3 * px_scale)]))

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
        fig.write_image(filename, scale=px_scale, width=width, height=height,
                        validate=True)

    return filename
