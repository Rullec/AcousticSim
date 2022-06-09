import sys
import os.path as osp

sys.path.append("../../dlls")
from softbody import softbody
import os

# os.chdir('D:\\Projects\\AcousticSim')
import scipy
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from matplotlib import pyplot as plt
from calc_stiff import calc_stiffness

SAMPLING_RATE = 44100 * 4


def get_sparse_K(K, dof):
    # val = [i for i in range(30)]
    # 1. rows
    # 2. cols
    # 3. data
    rows = K[0::3]
    cols = K[1::3]
    vals = K[2::3]
    # print(len(K) % 3)
    # print(len(rows))
    # print(len(cols))
    # print(len(vals))
    Kmat = coo_matrix((vals, (rows, cols)), shape=(dof, dof))
    # print(Kmat)
    return Kmat


def get_sparse_M(M_diag):
    dof = len(M_diag)
    rows = [i for i in range(dof)]
    cols = rows
    Mmat = coo_matrix((M_diag, (rows, cols)), shape=(dof, dof))
    return Mmat


def get_result(conf_path):

    print(f"[info] get result for {conf_path} cwd {os.getcwd()}")
    a = softbody()
    a.InitFromFile(conf_path)

    dampinga, dampingb = a.GetRayleighDamping()
    num_of_dof = len(a.GetVertexPos())

    # FEM stvk / neohookean elasticity
    # K = a.GetGlobalStiffnessMatrixEntrys()
    # for i in range(2, len(K), 3):
    #     K[i] *= -1
    # Kmat = get_sparse_K(K, num_of_dof)

    # linear elasticity
    Kmat = calc_stiffness(a)
    np.savetxt("log.py", Kmat.todense())
        

    # lumped diag M
    # M = a.GetMassMatrixDiag()
    # Mmat = get_sparse_M(M)
    # print(f"M {Mmat}")
    # raw symetric M
    M_entrys = a.GetGlobalRawMassMatrixEntrys()
    Mmat = get_sparse_K(M_entrys, num_of_dof)
    # print(f"M {Mmat}")

    # print(len(K))
    # print(len(M))

    # K_mat = coo_matrix(K,  (num_of_dof, num_of_dof))

    eigen_values, eigen_vecs = eigsh(A=Kmat, k=num_of_dof - 3, M=Mmat)
    eigen_values = np.real(eigen_values)
    eigen_vecs = np.real(eigen_vecs)
    # print("weight", weight)
    # plt.scatter(omega, weight)
    # plt.show()
    return eigen_values, eigen_vecs, dampinga, dampingb


def generate_sound(eigen_vals, eigen_vecs, dampa, dampb):
    # 1. get force
    print(f"dampa {dampa} dampb {dampb}")
    force = np.zeros(eigen_vecs.shape[0])
    force[0] = 0.5
    force[1] = 0.5
    force[2] = 0.5

    # 2. get modal equation RHS
    RHS = np.dot(eigen_vecs.T, force)
    interest_dof = 0

    # 3. get sound for each mode
    duration = 1.0
    t = np.linspace(0., duration, int(SAMPLING_RATE * duration))
    dt = 1.0 / (SAMPLING_RATE)
    sound = np.zeros(len(t))
    for _idx in range(len(eigen_vals)):
        '''
        qddot + (a + b k) qdot + k q = RHS
        k is eigen value

        solution:
        q(t) = dt * RHS / wd * e^{-1 * xi * w * t} * sin(wd * t)

        x(t) = eigen_vec * q
        '''
        eigen_val = eigen_vals[_idx]
        omega = np.sqrt(eigen_val)
        
        if omega > 20000 * 2 * np.pi or omega < 20 * 2 * np.pi or np.isnan(omega):
            continue
        else:
            # hearable!
            
            '''
            w = omega
            xi = (a + b * k) / (2 * w)
            wd = w * sqrt(1 - xi * xi)
            '''
            w = omega
            xi = (dampa + dampb * eigen_val) / (2 * w)
            wd = w * np.sqrt(1 - xi * xi)

            amp = RHS[_idx] / wd * np.exp(-1 * xi * w * t)
            wave = np.sin(wd * t) 

            print(f"mode {_idx} omega {omega:.1f} xi {xi:.1f} wd {wd:.1f} amp[0] {amp[0]:.1e}")            
            sound_cmp = wave * amp
            print(len(amp), len(wave), len(sound_cmp))
            sound += sound_cmp
    plt.plot(sound)
    plt.show()
    return sound


def save_to_wav(name, sound):
    from scipy.io.wavfile import write

    write(name, SAMPLING_RATE, 1e4 * sound)
    print(f"dump to {name}")
    # samplerate = 44100; fs = 100
    # plt.plot(sound)
    # plt.show()
    # t = np.linspace(0., 1., samplerate)
    # data = np.sin(2. * np.pi * fs * t)
    # write("example.wav", samplerate, data.astype(np.float))
    # print(f"dump sound to {name}")


if __name__ == "__main__":
    confs = [
        "config/obj_configs/obj0.json",
        # "config/obj_configs/obj38.json",
        # "config/obj_configs/obj98.json",
        # "config/obj_configs/beam/beam0.json",
        # "config/obj_configs/beam/beam1.json",
        # "config/obj_configs/beam/beam2.json",
        # "config/obj_configs/beam/beam3.json",
        # "config/obj_configs/beam/beam4.json",
        # "config/obj_configs/beam/beam5.json",
        # "config/obj_configs/beam/beam6.json",
        # "config/obj_configs/beam/beam7.json",
        # "config/obj_configs/beam/beam8.json",
        # "config/obj_configs/beam/beam9.json",
        # "config/obj_configs/beam/beam11.json",
        # "config/obj_configs/beam/beam12.json",
        # "config/obj_configs/beam/beam13.json",
        # "config/obj_configs/beam/beam14.json"
        # "config/obj_configs/plate/plate_122.json",
        # "config/obj_configs/plate/plate_434.json"
    ]

    for _idx, i in enumerate(confs):
        # plt.subplot(3, 4, _idx + 1)
        eigen_vals, eigen_vecs, dampa, dampb = get_result(i)
        name = osp.split(i)[-1]
        sound = generate_sound(eigen_vals, eigen_vecs, dampa, dampb)
        save_to_wav(name[:-5] + ".wav", sound)

        # for j in range(len(eigen_vals)):
        #     print(f"mode {j} omega {np.sqrt( eigen_vals[j])}")
        # plt.scatter(omega, weight, label=name, s=5)
        # plt.xlim([0, 1.2e5])
    # plt.legend()
    # plt.show()