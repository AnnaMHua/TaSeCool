import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


with open("Result_2020_08_23_18_08_1598222370.json") as f:
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
    graph._offsets3d = (x, y, z)
    title.set_text('gap={}'.format(gap))
    return graph,
x = [0]
y = [0]
z = [0]
# x, y, z = np.transpose(np.array(dic[0]))



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
title = ax.set_title('gap')
graph = ax.scatter(x, y, z, color='orange')
# x, y, z = np.transpose(np.array(dic[gapArray[0]]))
# graph._offsets3d = (x, y, z)

anim =FuncAnimation(fig, update_graph, gapArray,interval=1, blit=False )
plt.show()


