import json
import numpy as np
from viewer_3d import Viewer

convert_to_np_array_of_lst = lambda x: [np.array(i) for i in x]
air_density = 1.225


class cVibObj:
    def __init__(self):
        self.vertex_pos_lst = []
        self.vertex_normal_lst = []
        self.mode_info_lst = []  # 4 coef determine a mode
        self.vertex_coef_lst = []
        self.num_of_poles = 1
        self.pole_pos_lst = [np.zeros(3) for i in range(self.num_of_poles)]
        self.pole_weight_lst = [1 for i in range(self.num_of_poles)]

        # init calculation
        self.calc_vertex_sound_pressure()

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

    def calc_vertex_sound_pressure_of_each_mode(self):
        num_of_mode = self.mode_info_lst
        num_of_v = len( self.vertex_pos_lst)
        for mode_id_j in range(num_of_mode):
            '''
            in undamped case, mode vibration coef Aj = dt * (U^T finit)_j / wd_j
            '''
            for v_id_i in range(num_of_v):
                # calculate vertex sound pressure
                for local_dof in range(3):
                    # calculate vertex x,y,z accel
                    '''
                    \Delta a_{dof} = - U_{dof, j} * Aj * wd_j^2  * sin(wd_j * t)
                    '''
                    # ! need to get the whole eigen vecs ( for each dof )
                '''
                in undamped case, the sound pressure p_i(x) on vertex i for mode j is:
                \Delta a_i = - Uij * Aj * wd_j^2 * sin(wd_j * t)

                \Delta vec_a_i = [
                    \Delta a_(3 * i + 0),
                    \Delta a_(3 * i + 1),
                    \Delta a_(3 * i + 2),
                ]
                p_i(x) = rho * ni.dot(\Delta vec_a_i)
                '''
                mode_info = self.mode_info_lst[mode_id_j]
                coef_j, w, xi, wd = mode_info
                Aj = coef_j / wd



    def get_gradient(self):

        pass

    def get_energy(self):
        pass

    def optimize(self):

        pass


if __name__ == "__main__":
    path = "../../data/vib_steel_stvk.json"
    obj = cVibObj()
    obj.load(path)

    viewer = Viewer()
    print(obj.vertex_pos_lst)
    viewer.add_points(obj.vertex_pos_lst)
    num_of_v = len(obj.vertex_pos_lst)

    num_of_monopole = 1.0
    # 1. get gradient
    # 2. verify gradient
    # for i in range(num_of_v):
    #     st = obj.vertex_pos_lst[i]
    #     ed = obj.vertex_normal_lst[i] * 0.2
    #     viewer.add_arrow(st, ed)
    # viewer.show()