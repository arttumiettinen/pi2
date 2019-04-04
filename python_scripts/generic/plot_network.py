
import numpy as np
from mayavi import mlab


def draw_graph_simple(vertices, edges, plot_edges=True, point_size=0.03, edge_radius=0.01, color=(0.8, 0.0, 0.0), bgcolor=(1, 1, 1)):
    """
    Draws network graph, given vertex positions as Mx3 numpy array and edges as Nx2 numpy array that is an index list to vertices array.
    See example below.
    """

    mlab.figure(1, bgcolor=bgcolor)
    mlab.clf()

    pts = mlab.points3d(vertices[:, 0], vertices[:, 1], vertices[:, 2], scale_factor=point_size, scale_mode='none', color=color, resolution=5)

    if plot_edges:
        pts.mlab_source.dataset.lines = edges

        tube = mlab.pipeline.tube(pts, tube_radius=edge_radius)
        mlab.pipeline.surface(tube, color=color)

    mlab.show()




def example():
    """
    Example of simple network plotting functionality.
    """

    # Define vertex position array
    vertices = np.array([[0.0, 0.0, 0.0],
                         [1.0, 0.0, 0.0],
                         [0.0, 1.0, 0.0],
                         [0.0, 0.0, 1.0]])

    # Define edge array. Numbers are indices to vertices array.
    edges = np.array([[0, 1],
                      [0, 2],
                      [0, 3],
                      [1, 2],
                      [1, 3],
                      [2, 3]])

    # Draw
    draw_graph_simple(vertices, edges)


if __name__ == '__main__':
    example()
