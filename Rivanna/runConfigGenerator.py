'''

'''
import  json
import os

def LoadRunConfig(fname="runConfig.json"):
    data={}
    if os.path.isfile(fname):
        with open(fname) as fileIO:
            data = json.load(fileIO)
            # print(data)
            data['kzScan']['nbin']=100
    print(data)
    with open('test.json','w') as fileout:
        json.dump(data,fileout)

if __name__ == '__main__':
    LoadRunConfig()