from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Arrow3D(FancyArrowPatch):
    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)


def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    '''Add an 3d arrow to an `Axes3D` instance.'''

    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)


setattr(Axes3D, 'arrow3D', _arrow3D)


class Viewer():
    def __init__(self):
        fig = plt.figure()
        self.ax = Axes3D(fig)

    def add_points(self, pt_lst):
        x_lst = [i[0] for i in pt_lst]
        y_lst = [i[1] for i in pt_lst]
        z_lst = [i[2] for i in pt_lst]
        self.ax.scatter(x_lst, y_lst, z_lst)

    def add_arrow(self, st, ed):
        # ax.arrow3D(0,0,0,
        #    1,1,1,
        #    mutation_scale=20,
        #    arrowstyle="-|>",
        #    linestyle='dashed')
        self.ax.arrow3D(*st, *ed, mutation_scale=20, ec='green', fc='red')

    def show(self):
        plt.show()


if __name__ == "__main__":

    viewer = Viewer()
    viewer.add_arrow([0, 0, 0], [1, 1, 1])
    viewer.show()