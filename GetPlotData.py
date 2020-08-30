import glob
import json
import time

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

def ResultDumper(result=[],fname=""):
    if not fname:
        fname="./result/dHtri_PlotResult_d={}.json".format(d)

    with open(fname,"w") as fileio:
        json.dump(result,fileio)
    fileio.close()

if __name__ == '__main__':
    tic = time.time()
    d=10
    # print(selectWeyl(d))
    ResultDumper(selectWeyl(d))
    toc = time.time()
    print(toc-tic)