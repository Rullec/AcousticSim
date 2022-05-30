from matplotlib import pyplot as plt
import numpy as np
dt = 0 

def get_accel(disp, dt):
    x = []
    accel = []

    for i in range(1, len(disp) - 1):
        accel.append((disp[i - 1] + disp[i + 1] - 2 * disp[i]) / (dt * dt))
        x.append(i * dt)
    return x, accel


def get_val(file):
    global dt
    x_lst = []
    data_lst = []
    with open(file, 'r') as f:
        lines = f.readlines()
        dt = 1.0 / float(lines[0].split(" ")[0])
        lines = lines[1:]
        for _idx, line in enumerate(lines):
            # try:
            while line.find("  ") != -1:
                line = line.replace("  ", " ")
            data = [float(i.strip()) for i in line.split(" ")[5:]]

            x = [dt * i for i in range(len(data))]
            x_lst.append(x)
            data_lst.append(data)

            # plt.plot(x, data, label=f"{_idx}")
            # cycles = 0
            # for j in range(1, len(data) -1 ):
            #     if data[j] < data[j-1] and data[j] < data[j+1]:
            #         cycles+=1
            # print(f"cycles {cycles}")
        # break
        # except:
        #     pass
    return x_lst, data_lst
    # plt.legend()
    # plt.xlabel("time/s")
    # plt.ylabel("amp/m")
    # plt.show()


vx, vy = get_val("v_vib.txt")
mx, my = get_val("mode_vib.txt")
clamp = -1
plt.subplot(2, 2, 1)
for i in range(len(mx)):
    x = mx[i][:clamp]
    y = my[i][:clamp]
    plt.plot(x, y, label=f"mode {i}")
plt.legend()
plt.xlabel("time/s")

plt.subplot(2, 2, 2)
for i in range(len(vx)):
    x = vx[i][:clamp]
    y = vy[i][:clamp]
    plt.plot(x, y, label=f"v {i}")
plt.legend()
plt.xlabel("time/s")
plt.ylabel("amp/m")

plt.subplot(2, 2, 3)
for i in range(len(vx)):
    # x = vx[i][:clamp]

    x, v_accel = get_accel(vy[i][:clamp], dt)
    plt.plot(x, v_accel, label=f"v_acc {i}")
plt.legend()
plt.xlabel("time/s")
plt.ylabel("accel [m.s-2]")


plt.subplot(2, 2, 4)

play_x, play_y = get_val("wave_play.txt")

plt.plot(play_x[0], play_y[0], label=f"play")
plt.legend()
plt.xlabel("time/s")
plt.ylabel("pressure pa")

plt.show()
