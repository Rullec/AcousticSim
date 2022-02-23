from matplotlib import pyplot as plt
import numpy as np

with open("vib.txt", 'r') as f:
    for _idx, line in enumerate(f.readlines()):
        # try:
            while line.find("  ") != -1:
                line = line.replace("  ", " ")
            data = [float(i.strip()) for i in line.split(" ")[5:]]

            dt = 2e-3
            x = [dt * i for i in range(len(data))]
            plt.plot(x, data, label=f"{_idx}")
            cycles = 0
            for j in range(1, len(data) -1 ):
                if data[j] < data[j-1] and data[j] < data[j+1]:
                    cycles+=1
            print(f"cycles {cycles}")
            # break
        # except:
        #     pass
    plt.legend()
    plt.xlabel("time/s")
    plt.ylabel("amp/m")
    plt.show()