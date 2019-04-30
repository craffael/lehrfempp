from sys import argv

import matplotlib.pyplot as plt
import numpy as np

if len(argv) < 2:
    print('usage: python plot_mesh.py writeMatplotlib_output.txt')
    exit(-1)

points = dict()
segments = dict()
triangles = dict()
quadrilaterals = dict()

with open(argv[1]) as f:
    for line in f:
        data = line.rstrip().split(' ')
        geometry_type = data[0]
        index = int(data[1])
        coordinates = [
            np.array([float(x), float(y)])
            for x, y in zip(data[2::2], data[3::2])
        ]

        if geometry_type == 'Point':
            geometry_list = points
        elif geometry_type.startswith('Segment'):
            geometry_list = segments
        elif geometry_type.startswith('Tria'):
            geometry_list = triangles
        elif geometry_type.startswith('Quad'):
            geometry_list = quadrilaterals
        else:
            raise RuntimeError('Unknown geometry')

        geometry_list[index] = coordinates

x_min, x_max = float('inf'), float('-inf')
y_min, y_max = float('inf'), float('-inf')

# plot vertices
for idx, coords in points.items():
    for coord in coords:
        x_min = min(coord[0], x_min)
        x_max = max(coord[0], x_max)
        y_min = min(coord[1], y_min)
        y_max = max(coord[1], y_max)

    plt.annotate(
        idx,
        coords[0],
        fontsize='xx-large',
        weight='bold',
        bbox=dict(boxstyle='circle', facecolor='none', edgecolor='red'),
        color='r',
        ha='center',
        va='center'
    )


def poly_coeffs(x, coeffs):
    order = len(coeffs)
    y = 0

    for i in range(order):
        y += coeffs[i] * x ** (order - (i + 1))

    return y


def collinear(p0, p1, p2):
    x1, y1 = p1[0] - p0[0], p1[1] - p0[1]
    x2, y2 = p2[0] - p0[0], p2[1] - p0[1]
    return abs(x1 * y2 - x2 * y1) < 1e-12


# plot cells
for idx, coords in segments.items():
    for coord in coords:
        x_min = min(coord[0], x_min)
        x_max = max(coord[0], x_max)
        y_min = min(coord[1], y_min)
        y_max = max(coord[1], y_max)

    # SegmentO1
    if len(coords) == 2:
        coordinates = np.vstack(coords)
        x_coords = coordinates[:, 0]
        y_coords = coordinates[:, 1]
        plt.plot(x_coords, y_coords, 'k-')

        plt.annotate(
            idx,
            [x_coords.mean(), y_coords.mean()],
            fontsize='xx-large',
            ha='center',
            va='center'
        )

    # SegmentO2
    elif len(coords) == 3:
        vertex_0 = coords[0]
        vertex_1 = coords[1]
        midpoint = coords[2]

        # check if points are collinear
        if collinear(vertex_0, midpoint, vertex_1):
            plt.plot(
                [vertex_0[0], vertex_1[0]], [vertex_0[1], vertex_1[1]], 'k-'
            )
        else:
            coefficients = np.polyfit(
                [vertex_0[0], midpoint[0], vertex_1[0]],
                [vertex_0[1], midpoint[1], vertex_1[1]],
                2
            )
            x = np.linspace(vertex_0[0], vertex_1[0], 1000)
            plt.plot(x, poly_coeffs(x, coefficients), 'k-')

        plt.annotate(
            idx,
            midpoint,
            fontsize='xx-large',
            ha='center',
            va='center'
        )
# plot cells
cell_data = {'limegreen': triangles, 'm': quadrilaterals}

for color, cells in cell_data.items():
    for idx, coords in cells.items():
        for coord in coords:
            x_min = min(coord[0], x_min)
            x_max = max(coord[0], x_max)
            y_min = min(coord[1], y_min)
            y_max = max(coord[1], y_max)

        plt.annotate(
            idx,
            np.vstack(coords).mean(0),
            fontsize='xx-large',
            color=color,
            weight='bold',
            ha='center',
            va='center'
        )

offset = 1e-2
plt.xlim(left=x_min - offset, right=x_max + offset)
plt.ylim(bottom=y_min - offset, top=y_max + offset)
plt.axis('off')
plt.show()
