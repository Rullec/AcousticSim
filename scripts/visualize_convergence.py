from matplotlib import pyplot as plt
name_lst = ["../log.150", "../log.500", "../log.1000"]
for  _idx, name in enumerate(name_lst):
    with open(name) as f:
        lines = [line for line in f.readlines() if line.find("res max") != -1]
        value_lst = []
        for line in lines:
            try:
                value_lst.append(float(line.split(" ")[-1]))
            except:
                pass
        plt.subplot(1, 3, _idx + 1)
        plt.plot(value_lst)
        plt.title(f"dim = %s" % name[name.find("log.") + 4:])
        plt.ylim([0, 10])
plt.show()