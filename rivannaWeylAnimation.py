import glob
import json

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib.animation import FuncAnimation
from math import pi

def selectWeyl(d):
    path = "/Users/mhua/Documents/Research/TaSe/myNotes/TaSeCool/result/dHtri_Result_*.json"
    dic = {}
    for filename in glob.glob(path):
        with open(filename,'r') as f:
            data = json.load(f)
            for value in data:
                if value["HamiVal"][1] - value["HamiVal"][0] < d:
                    if value["gap"] in dic.keys():
                        dic[value["gap"]].append(value["Pos"])
                    else:
                        dic[value["gap"]] = [value["Pos"]]

            f.close()


    return dic

def update_graph(gap):
    x, y, z = np.transpose(np.array(dic[gap]))
    graph.set_data(x, y)
    graph.set_3d_properties(z)

    title.set_text('3D Test, time={}'.format(gap))
    return title, graph,

if __name__ == '__main__':
    d=10
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    title = ax.set_title('3D Test')

    plt.xlim(0, pi)
    plt.ylim(0, pi)
    ax.set_zlim3d(0,pi)

    # initialize scatters
    dic = selectWeyl(d)
    print(list(dic.keys()))
    x, y, z = np.transpose(np.array(dic[list(dic.keys())[0]]))
    graph, = ax.plot(x, y, z, linestyle="", marker="o")

    anim = FuncAnimation(fig, update_graph, list(dic.keys()), interval=1, blit=False)
    anim.save("dHriWeyl3Danim_d=" + str(d) + ".mp4")


