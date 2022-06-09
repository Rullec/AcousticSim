import sys
import numpy as np
import os.path as osp

sys.path.append("../../dlls")
from softbody import softbody
import os


def calc_sign_volume(tet_vertex_pos_lst):
    M = np.zeros([4, 4])
    M[:, 0] = 1
    for i in range(4):
        # print(f"v{i} pos {tet_vertex_pos_lst[i]}")
        M[i, 1:] = tet_vertex_pos_lst[i]
    Mdet_6V = np.linalg.det(M)
    return Mdet_6V / 6


os.chdir('../../')


def calcB(tet_vertex_pos_lst):
    # 1. calculate coef
    M = np.zeros([4, 4])
    M[:, 0] = 1
    for i in range(4):
        M[i, 1:] = tet_vertex_pos_lst[i]
    Mdet_6V = np.linalg.det(M)
    assert Mdet_6V > 0
    C = Mdet_6V * np.linalg.inv(M).T

    beta = C[:, 1]
    gamma = C[:, 2]
    delta = C[:, 3]

    # 2. calculate [B]
    B = np.zeros([6, 12])
    for i in range(4):
        B[0, 3 * i] = beta[i]
        B[1, 3 * i + 1] = gamma[i]
        B[2, 3 * i + 2] = delta[i]
    for i in range(4):
        # row 3
        B[3, 3 * i] = gamma[i]
        B[3, 3 * i + 1] = beta[i]

        # row 4
        B[4, 3 * i + 1] = delta[i]
        B[4, 3 * i + 2] = gamma[i]

        # row 5
        B[5, 3 * i + 0] = delta[i]
        B[5, 3 * i + 2] = beta[i]
    B = B / Mdet_6V
    # B = np.round(B, 1)
    # print(f"B \n{B}")
    # exit()
    return B


def calcD(body):
    nu = body.GetPoissonRatio()
    E = body.GetYoungsModulus()
    print(f"poisson ratio {nu}")
    print(f"E {E}")

    D = np.zeros([6, 6])
    D[:3, :3] = np.array([
        [1 - nu, nu, nu],
        [nu, 1 - nu, nu],
        [nu, nu, 1 - nu],
    ])
    for j in range(3, 6):
        D[j, j] = (1 - 2 * nu) / 2
    # print(f"D first {D}")
    D *= E / (1 + nu) / (1 - 2 * nu)
    # print(f"D final {D}")
    # exit()
    return D


from scipy.sparse import coo_matrix

from tqdm import tqdm


def calc_stiffness(body: softbody):

    # 1. calculate D
    D = calcD(body)

    vertex_lst = body.GetVertexPos()

    tet_vid_lst = body.GetTetVertexIdx()
    row_lst = []
    col_lst = []
    val_lst = []
    num_of_tet = len(tet_vid_lst)
    print(f"[info] total num of tet {num_of_tet}")
    dof = len(vertex_lst)

    for t in tqdm(range(num_of_tet)):
        # 2. calculate B for each element
        tet_vpos = [vertex_lst[3 * vid:3 * vid + 3] for vid in tet_vid_lst[t]]
        B = calcB(tet_vpos)
        # print(tet_vpos)
        # exit()
        V_sign = calc_sign_volume(tet_vpos)
        V_abs = np.abs(V_sign)
        # print(f"V sign and abs {V_sign} {V_abs}")
        # 3. generate K
        # K = V * BT * D * B
        K_ele = V_abs * np.dot(np.dot(B.T, D), B)
        for a in range(4):
            va = tet_vid_lst[t][a]
            for b in range(4):
                vb = tet_vid_lst[t][b]
                for j in range(3):
                    for k in range(3):
                        global_row_id = 3 * va + j
                        global_col_id = 3 * vb + k
                        row_lst.append(global_row_id)
                        col_lst.append(global_col_id)
                        val_lst.append(K_ele[3 * a + j, 3 * b + k])

    K_global = coo_matrix((val_lst, (row_lst, col_lst)), shape=(dof, dof))
    return K_global


if __name__ == "__main__":
    np.set_printoptions(suppress=True)
    conf_path = "config/obj_configs/obj0.json"
    body = softbody()
    body.InitFromFile(conf_path)
    calc_stiffness(body)