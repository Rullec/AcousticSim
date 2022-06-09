import json
import numpy as np
from viewer_3d import Viewer

convert_to_np_array_of_lst = lambda x: [np.array(i) for i in x]


class cVibObj:
    def __init__(self):
        self.vertex_pos_lst = []
        self.vertex_normal_lst = []
        self.mode_info_lst = []  # 4 coef determine a mode
        self.vertex_coef_lst = []

    def load(self, path):
        with open(path) as f:
            cont = json.load(f)

        coef = cont["coef"]
        geo = cont["geometry"]
        modes = cont["modes"]

        # 1. parse coef
        self.vertex_coef_lst = convert_to_np_array_of_lst(coef)
        print(len(self.vertex_coef_lst))
        # 2. parse geo
        self.vertex_normal_lst = convert_to_np_array_of_lst(geo["normal_lst"])
        self.vertex_pos_lst = convert_to_np_array_of_lst(geo["pos_lst"])
        print(len(self.vertex_pos_lst))

        # 3. parse modes
        self.mode_info_lst = convert_to_np_array_of_lst(modes)
        print(len(self.mode_info_lst))


if __name__ == "__main__":
    path = "vib_steel_stvk.json"
    obj = cVibObj()
    obj.load(path)

    viewer = Viewer()
    print(obj.vertex_pos_lst)
    viewer.add_points(obj.vertex_pos_lst)
    num_of_v = len(obj.vertex_pos_lst)

    for i in range(num_of_v):
        st = obj.vertex_pos_lst[i]
        ed = obj.vertex_normal_lst[i] * 0.2
        viewer.add_arrow(st, ed)
    viewer.show()