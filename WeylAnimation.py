import json
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib.animation import FuncAnimation

x = [0]
y = [0]
z = [0]


with open("./result/Result_2020_08_26_14_08_1598468252.json") as f:
    data = json.load(f)
    dic = {}
    gapArray=[]
    for key, value in data.items():
        dic[value["gap"]] = []
    for key, value in data.items():
        dic[value["gap"]].append(value["Pos"])
    gapArray = list(dic.keys())


def update_graph(gap):
    x, y, z = np.transpose(np.array(dic[gap]))
    graph.set_data(x, y)
    graph.set_3d_properties(z)

    title.set_text('3D Test, time={}'.format(gap))
    return title, graph,


fig=plt.figure()
ax=p3.Axes3D(fig)
title = ax.set_title('3D Test')


# initialize scatters
x, y, z = np.transpose(np.array(dic[gapArray[0]]))
graph, = ax.plot(x, y, z, linestyle="", marker="o")


anim =FuncAnimation(fig, update_graph, gapArray,interval=1, blit=False )
plt.show()
