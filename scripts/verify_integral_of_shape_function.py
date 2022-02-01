from cgi import test
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
    
def judge_point_inside_tet(tet_pt_lst, cur_pt):
    # https://blog.csdn.net/ycl295644/article/details/49615377
    D0 = np.array(
        [
            [*(tet_pt_lst[0]), 1],
            [*(tet_pt_lst[1]), 1],
            [*(tet_pt_lst[2]), 1],
            [*(tet_pt_lst[3]), 1]
        ], dtype = np.float32
    )
    sign_det_D0 = np.sign(np.linalg.det(D0))
    # print(f"----begin to judge---")
    # print(f"sign_det_D0 {sign_det_D0}")
    for j in range(4):
        Dj = deepcopy(D0)
        Dj[j, :3] = cur_pt

        sign_det_Dj = np.sign(np.linalg.det(Dj))
        if sign_det_Dj == 0 :
            continue
        if sign_det_Dj != sign_det_D0:
            # outside of this tet
            # print(f"sign_det_D{j} {sign_det_Dj}")
            # print(f"----end to judge---")
            return False
    # print(f"----end to judge---")
    return True

def get_box(A, B, C, D):
    pt_collection = np.array([A, B, C, D])
    min_x, min_y, min_z = np.min(pt_collection, axis = 0)
    max_x, max_y, max_z = np.max(pt_collection, axis = 0)
    # print(min_x, min_y, min_z)
    # print(max_x, max_y, max_z)
    return min_x, max_x - min_x, min_y, max_y - min_y, min_z, max_z - min_z

def test_judge_point_inside_tet(tet_pt_lst):
    _tet_pt_lst = np.array(tet_pt_lst)
    min_x, x_range, min_y, y_range, min_z, z_range  = get_box(_tet_pt_lst[0], _tet_pt_lst[1], _tet_pt_lst[2], _tet_pt_lst[3])
    # = get_box(A, B, C, D)
    # random points in a unit cube
    pts = np.random.rand(10000, 3)
    pts[:, 0] = pts[:, 0] * x_range + min_x
    pts[:, 1] = pts[:, 1] * y_range + min_y
    pts[:, 2] = pts[:, 2] * z_range + min_z
    fig = plt.figure()
    ax = Axes3D(fig)
    new_pts = []
    for j in range( pts.shape[0]):
        cur_pt = pts[j]
        if judge_point_inside_tet(tet_pt_lst, cur_pt) == True:
            new_pts.append(cur_pt)
        else:
            print(f"ignore")
    new_pts = np.array(new_pts)
    print(new_pts.shape)
    ax.scatter(new_pts[:, 0], new_pts[:, 1], new_pts[:, 2])
    
    plt.show()
    
    
class cShapeFunc:
    def __init__(self, a, b, c, d, pt_lst) -> None:

        
        self.a = a
        self.bcd = np.array([b, c, d])
        self.pt_lst = pt_lst
    
    def evaluation(self, pt):
        assert len(pt) == 3
        if judge_point_inside_tet(self.pt_lst, pt) == True:
            # inside tet: shape function valid
            
            return self.a + np.dot(self.bcd, pt)
        else:
            # outside the tet
            # print(f"judge cur pt {pt} is outside of tet {self.pt_lst}")
            return 0
    

        
class cTetShapeFunc:
    def __init__(self, A, B, C, D) -> None:
        self.pt_lst = [A, B, C, D]
        self.__init_shape_function__()
        
    def __init_shape_function__(self):
        D = np.ones([4, 4])
        for j in range(4):
            D[j, 1:4] = self.pt_lst[j]
        det_D = np.linalg.det(D)
        # print(f"get D = {D}, det(D) = {det_D}")
        if np.abs(det_D) < 1e-6:
            raise ValueError("det(D)  abs < 1e-6, ill-conditioned!")
        
        self.shape_func_lst = []
        for i in range(4):
            # print(f"begin to build N{i}")
            '''
            given N_i(x, y, z) = a_i + b_i * x + c_i * y + d_i * z
            solve for a_i, b_i, c_i, d_i
            a_i = det(D_0) / det(D) ; 其中D_0是D的第0列被替换为b向量的矩阵
            b_i = det(D_1) / det(D) ； 其中D_1是D的第1列被替换为b向量的矩阵
            c_i = det(D_2) / det(D) ； 其中D_2是D的第2列被替换为b向量的矩阵
            d_i = det(D_3) / det(D) ； 其中D_3是D的第3列被替换为b向量的矩阵
            '''
            b = np.zeros(4)
            b[i] = 1
            # print(f"for P{i}, b = {b}")
            # build ai, bi, ci, di; 
            coef_lst = []
            for j in range(4):
                # calculate Dj
                D_j = deepcopy(D)
                D_j[:, j] = b
                # print(f"D{j} = {D_j}")
                det_Dj = np.linalg.det(D_j)
                coef_lst.append(det_Dj / det_D)
            ai, bi, ci, di = coef_lst
            # print(f"ai = {ai}, bi = {bi}, ci = {ci}, di = {di}")
            Ni = cShapeFunc(ai, bi, ci, di, self.pt_lst)
            self.shape_func_lst.append(Ni)
            # exit()
            

            
            

def calc_tet_volume(A, B, C, D):
    '''
    calculate tetrahedraon volume
    '''
    AB = B - A
    AC = C - A
    AD = D - A
    V =  1.0 / 6 * (np.cross(AB, AC)).dot(AD)
    return V

def calc_J(A, B, C, D):
    J = np.zeros([3, 3])
    
    J[:, 0] = B - A 
    J[:, 1] = C - A 
    J[:, 2] = D - A 
    return J


def calc_value_for_ninj(shape_func_lst, i, j, cur_pt):
    Ni = shape_func_lst[i]
    Nj = shape_func_lst[j]
    return Ni.evaluation(cur_pt) * Nj.evaluation(cur_pt)
from tqdm import tqdm

def mc_est(A, B, C, D, shape_func_lst):
    num_of_pts = 10000
    min_x, x_range, min_y, y_range, min_z, z_range = get_box(A, B, C, D)
    # random points in a unit cube
    pts = np.random.rand(num_of_pts, 3)
    pts[:, 0] = pts[:, 0] * x_range + min_x
    pts[:, 1] = pts[:, 1] * y_range + min_y
    pts[:, 2] = pts[:, 2] * z_range + min_z
    prob = 1.0 / (x_range * y_range * z_range)
    data =  []

    for i in range(4):
        for j in tqdm(range(4)):
            total_integral = 0        

            for _idx in range(num_of_pts):
                cur_pt = pts[_idx]
                
                value = calc_value_for_ninj(shape_func_lst, i, j, cur_pt)
                total_integral += value / prob 
            
            total_integral /= num_of_pts
            data.append((i, j, total_integral))
    return data

def test_shape_func_01(shape_func_lst, pt_lst):
    assert len(pt_lst) == len(shape_func_lst)
    assert len(pt_lst) == 4
    for i in range(4):
        Ni = shape_func_lst[i]
        for pt_idx in range(4):
            cur_pt = pt_lst[pt_idx]
            # print(f"cur pt {cur_pt}")
            eval = Ni.evaluation(cur_pt)
            assert eval == (i == pt_idx), f"eval {eval} for N{i} in idx{pt_idx}"

if __name__ == "__main__":
    A = np.array([0, 0, 0])
    B = np.array([1, 0, 0])
    C = np.array([0, 1, 0])
    D = np.array([0, 0, 6])
    pt_lst=  [A, B, C, D]
    # test_judge_point_inside_tet([A, B, C, D])
    
    tet = cTetShapeFunc(A, B, C, D)
    shape_func_lst = tet.shape_func_lst
    test_shape_func_01(shape_func_lst, pt_lst)
    
    # exit()

    # volume is 1
    volume = calc_tet_volume(A, B, C, D)
    # print(f"volume {volume}")
    # 1. get the analytic result
    J = calc_J(A, B, C, D)
    det_J = np.linalg.det(J)
    # print(f"J {J} \ndet {det_J}")
    integral_same = det_J / 60
    integral_diff = det_J / 120
    print(f"integral same {integral_same} integral diff {integral_diff}")

    # 2. monte carlo estimation
    data_lst = mc_est(A, B, C, D, shape_func_lst)
    for data in data_lst:
        i, j, val = data
        print(f"MonteCarlo result: \int N{i}*N{j} = {val}")