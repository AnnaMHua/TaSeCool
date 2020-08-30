import TaSeHamiltonian as Tase
from numpy import linalg as la
import os
import sys
import json
import math
import numpy as np
from multiprocessing.pool import ThreadPool
import multiprocessing
import  logging
import datetime
from math import cos, sin, sqrt, pi
import time


class WeylFinder(object):
    # Scan range
    krScanMin = 0
    krScanNBin =10
    krScanMax = pi

    ksScanMin = 0
    ksScanNBin = 10
    ksScanMax = pi
    ksStep = pi/ksScanNBin


    kzScanMin = 0
    kzScanNBin = 10
    kzScanMax = pi
    kzStep = pi/kzScanNBin


    gapScanMin = 0
    gapScanNBin = 50
    gapScanMax = 200

    jobnumber = 0


    def __init__(self, HamiClass):
        if os.path.isfile("./resource/deve.infor"):
            with open("./resource/deve.infor",encoding="utf8") as infor:
                chatInfor=infor.read()
            print(chatInfor)
        self.HamiClass = HamiClass


    def jobArrayto3D(self,jobnum):
        """
        transform job array number to kr,ks,kz idex
        """

        # number of cuts in each direction
        cutN = 10
        if 0<=jobnum <1000:
            self.jobnumber = jobnum
        else:
            raise IOError

        krStep = pi / cutN

        rindex = jobnum // 100
        sindex = (jobnum-100*rindex) //10
        zindex = jobnum -100*rindex -10*sindex
        # print(rindex, sindex, zindex)

        self.krScanMin = rindex * krStep
        self.krScanMax = (rindex + 1) * krStep

        self.ksScanMin = sindex * krStep
        self.ksScanMax = (sindex + 1) * krStep


        self.kzScanMin = zindex * krStep
        self.kzScanMax = (zindex + 1) * krStep
        # print(self.krScanMin,self.krScanMax)


    def HamiSingleScaner(self, SpacialPosgap=[]):
        # [a,b,c] = [(SpacialPos[0] -SpacialPos[1])/sqrt(2),(SpacialPos[0] +SpacialPos[1])/sqrt(2),SpacialPos[2]]
        # if self.HamiClass.dim ==16:
        val = la.eigvalsh(self.HamiClass.HamiltonianMatrix(spacialPos=SpacialPosgap[0:3],gap=SpacialPosgap[-1]))

        # return list(SpacialPosgap)
        return {"Pos":SpacialPosgap[0:3].tolist(),"gap":SpacialPosgap[-1],"HamiVal":val.tolist()}

    def HScanner(self):
        krRange = np.linspace(self.krScanMin, self.krScanMax, self.krScanNBin)
        ksRange = np.linspace(self.ksScanMin, self.ksScanMax, self.ksScanNBin)
        kzRange = np.linspace(self.kzScanMin, self.kzScanMax, self.kzScanNBin)
        gapRange = np.linspace(self.gapScanMin, self.gapScanMax, self.gapScanNBin)
        scanPosGap = np.reshape(np.transpose(np.meshgrid(krRange,ksRange,kzRange,gapRange)),(-1,4))
        scanResult = list(map(self.HamiSingleScaner,scanPosGap))
        # print(len(scanResult))
        # print(scanResult)
        self.ResultDumper(result=scanResult)
        # return  scanResult
        return None

    def ResultDumper(self,result=[],fname=""):
        '''

        '''
        if not fname:
            fname="./result/dHtri_Result_{}.json".format(self.jobnumber)

        with open(fname,"w") as fileio:
            json.dump(result,fileio)
        fileio.close()




if __name__ == '__main__':
    tic = time.time()
    test = WeylFinder(Tase.Henlarge("trivial"))
    print(test.HamiClass.dim)
    test.jobArrayto3D(int(sys.argv[1]))
    test.HScanner()
    toc = time.time()
    print(toc-tic)

