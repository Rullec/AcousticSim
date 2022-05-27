import numpy as np
from matplotlib import pyplot as plt

np.random.seed(0)

np.set_printoptions(suppress=True)
u0 = np.array([0.0, 0.0])
u1 = np.array([1.0, 0.5])
u2 = np.array([0.3, 2.0])

x0 = np.array([0.0, 0.0])
x1 = np.array([3.0, 1])
x2 = np.array([0.9, 4.0])


def normalize(vec):
    return np.array(vec) / np.linalg.norm(vec)


def calcDs(x0, x1, x2):
    dx0 = x1 - x0
    dx1 = x2 - x0
    dx = np.zeros([2, 2])
    dx[:, 0] = dx0
    dx[:, 1] = dx1
    return dx


def calcDm(u0, u1, u2):
    du0 = u1 - u0
    du1 = u2 - u0
    du = np.zeros([2, 2])
    du[:, 0] = du0
    du[:, 1] = du1
    return du


def CalcF(u0, u1, u2, x0, x1, x2):
    Ds = calcDs(x0, x1, x2)
    Dm = calcDm(u0, u1, u2)
    F = np.matmul(Ds, np.linalg.inv(Dm))
    return F


def CalcTriS(u0, u1, u2):
    a = u1 - u0
    b = u2 - u0

    return 0.5 * np.abs(np.cross(a, b))


F = CalcF(u0, u1, u2, x0, x1, x2)
Dm = calcDm(u0, u1, u2)
Dminv = np.linalg.inv(Dm)

print(f"Dm {Dm}")
Dm_det = np.linalg.det(Dm)
tri_area = CalcTriS(u0, u1, u2)
print(f"Dm_det {Dm_det}")
print(f"tri area {tri_area}")
print(f"Dminv {Dminv}")

exit()


def calcEnergy(u0, u1, u2, x0, x1, x2):
    F = CalcF(u0, u1, u2, x0, x1, x2)
    F0 = F[:, 0]
    F1 = F[:, 1]
    F0_norm = np.linalg.norm(F0)
    F1_norm = np.linalg.norm(F1)
    C = [F0_norm - 1, F1_norm - 1]
    E = np.dot(C, C)
    return E


def draw_tri(x0, x1, x2, name):
    x = [x0[0], x1[0], x2[0], x0[0]]
    y = [x0[1], x1[1], x2[1], x0[1]]
    plt.plot(x, y, label=name)


# F = CalcF(u0, u1, u2, x0, x1, x2)
# print(f"F {F}, {np.linalg.norm( F[:,0])}")
def rotmat2d(theta):
    m = np.identity(2)
    cosz = np.cos(theta)
    sinz = np.sin(theta)
    m[0, 0] = cosz
    m[0, 1] = -sinz
    m[1, 0] = sinz
    m[1, 1] = cosz
    return m


def calcpAndq(u0, u1, u2, x0, x1, x2):
    # 1. calc F
    F = CalcF(u0, u1, u2, x0, x1, x2)
    F0 = F[:, 0]
    F1 = F[:, 1]
    F0_norm = np.linalg.norm(F0)
    F1_norm = np.linalg.norm(F1)

    # 2. calc a, b, c, d
    Dminv = np.linalg.inv(calcDm(u0, u1, u2))
    a = Dminv[0, 0]
    c = Dminv[0, 1]
    b = Dminv[1, 0]
    d = Dminv[1, 1]

    # 3. calc vec m = \Delta x0 and n = \Delta x1
    Ds = calcDs(x0, x1, x2)
    m = Ds[:, 0]
    n = Ds[:, 1]

    # 4. calculate p and q
    '''
    part1 = (|F0| - 1) / |F0| * (a * m + b * n)
    part2 = (|F1| - 1) / |F1| * (c * m + d * n)
    p = part1 * a + part2 * c
    q = part1 * b + part2 * d
    '''
    part1 = (F0_norm - 1) / F0_norm * (a * m + b * n)
    part2 = (F1_norm - 1) / F1_norm * (c * m + d * n)
    p = part1 * a + part2 * c
    q = part1 * b + part2 * d
    return p, q


cur_e = calcEnergy(u0, u1, u2, x0, x1, x2)
print(f"cur e = {cur_e}")


def calcIncre(p, q, length=1e-3):
    dx0 = length * normalize([1, -1])
    '''
    p.dx0 + q.dx1 = 0
    let cos \theta = [-1, 1]
    |q||dx1| * cos theta = -p.dx0

    |dx1| = -p.dx0 / (cos_theta * |q|)
    '''
    q_length = np.linalg.norm(q)
    # cos_theta = (np.random.rand() - 0.5) * 2
    cos_theta = 0.99
    # print(f"cos theta = {cos_theta}")
    dx1_length = -np.dot(p, dx0) / (cos_theta * q_length)
    # print(f"dx1_length {dx1_length}")

    q_normalized = normalize(q)
    theta = np.arccos(cos_theta)
    # assert theta == 0
    # print(f"q_normalized {q_normalized}")
    dx1 = np.matmul(rotmat2d(theta), q_normalized) * dx1_length

    verify = np.dot(p, dx0) + np.dot(q, dx1)
    # print(f"dx1 {dx1}")
    # print(f"verify {verify}")
    return dx0, dx1


# dt = 1e-3
# x1_lst = [[], []]
# x2_lst = [[], []]
# tri_lst = []
# for _ in range(500):
#     tri_lst.append([x0, x1, x2])

#     p, q = calcpAndq(u0, u1, u2, x0, x1, x2)

#     dx0, dx1 = calcIncre(p, q, dt)
#     # print(f"dx0 {dx0}, dx1 {dx1}")
#     # dx0 = np.ones(2) * dt
#     # dx1 = np.ones(2) * dt
#     x1 = dx0 + x1
#     x2 = dx1 + x2
#     new_E = calcEnergy(u0, u1, u2, x0, x1, x2)
#     print(f"new_E {new_E}")

#     x1_lst[0].append(x1[0])
#     x1_lst[1].append(x1[1])

#     x2_lst[0].append(x2[0])
#     x2_lst[1].append(x2[1])

# # plt.plot(x1_lst[0], x1_lst[1], label = "x1")
# # plt.plot(x2_lst[0], x2_lst[1], label = "x2")
# gap = 10
# for i in range(gap):
#     draw_tri(* (tri_lst[len(tri_lst) // gap * i]), f"{i}-th")

# plt.legend()
# plt.show()