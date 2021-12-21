import pyvista as pv
import tetgen
import numpy as np

pv.set_plot_theme('document')

sphere = pv.Sphere()
tet = tetgen.TetGen(sphere)
# tet.make_manifold()
tet.tetrahedralize(order=1, mindihedral=20, minratio=1.5)
grid = tet.grid
# grid.save("save.cdb")
# grid.plot(show_edges=True)
# print(grid)
cells = grid.cells


def parse_cells(cells_array):
    idx = 0
    cell_node_lst = []
    while idx < len(cells_array):
        num_of_pts_in_this_node = cells_array[idx]
        print(num_of_pts_in_this_node)
        pt_lst = [
            cells_array[idx + j + 1] for j in range(num_of_pts_in_this_node)
        ]
        cell_node_lst.append(pt_lst)
        print(f"cur pt lst {pt_lst}")
        idx += num_of_pts_in_this_node + 1
        exit()


parse_cells(cells)
print(f"cells {type(cells)}")
print(f"len cells {len(cells)}")
print(f"len cells {cells.shape}")
exit()
# for cell in:
#     print(cell)
print(f"cells {grid.cells}")

print(f"points {grid.points}")
exit()

# get cell centroids
cells = grid.cells.reshape(-1, 5)[:, 1:]
cell_center = grid.points[cells].mean(1)

# extract cells below the 0 xy plane
mask = cell_center[:, 2] < 0
cell_ind = mask.nonzero()[0]
subgrid = grid.extract_cells(cell_ind)

# advanced plotting
plotter = pv.Plotter()
plotter.add_mesh(subgrid, 'lightgrey', lighting=True, show_edges=True)
plotter.add_mesh(sphere, 'r', 'wireframe')
plotter.add_legend([[' Input Mesh ', 'r'], [' Tessellated Mesh ', 'black']])
plotter.show()

# Compute cell quality0

cell_qual = subgrid.compute_cell_quality()['CellQuality']

# Plot quality

subgrid.plot(scalars=cell_qual,
             stitle='Quality',
             cmap='bwr',
             clim=[0, 1],
             flip_scalars=True,
             show_edges=True)
