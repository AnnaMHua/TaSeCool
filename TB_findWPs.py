'''
Multi thread Verison Wire Point Scanner

Autor: Siyu Jian

'''
from typing import IO

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

from tqdm import *

logger=logging.getLogger(__name__)
logger.setLevel(0)

class WireFinder(object):
    '''

    '''

    #Multi-Thread setting
    MaximumThread=8
    MaximumJobs=8    # how many sub jobs
    useSysCores=True  # use system detected number of cores instead of user specified

    #Scan range
    krScanMin=0
    krScanNBin=100
    krScanMax = pi/sqrt(2)

    ksScanMin = 0
    ksScanNBin = 100
    ksScanMax = pi/sqrt(2)

    kzScanMin = 0
    kzScanNBin = 100
    kzScanMax = pi

    gapScanMin = 0
    gapScanNBin = 50
    gapScanMax = 200

    saveFname=""
    _finishedWorkCount=0
    _finishedWorkCount_Update = 0

    _progressBarArray={}


    def __init__(self):
        if os.path.isfile("./resource/deve.infor"):
            with open("./resource/deve.infor",encoding="utf8") as infor:
                chatInfor=infor.read()
            print(chatInfor)

    def LoadRunConfig(self,fname="runConfig.json"):
        if os.path.isfile(fname):
            with open(fname) as fileIO:
                data=json.load(fileIO)
                try:
                    self.useSysCores=bool(data["runConfig"]["useSysCores"])
                    # print(self.useSysCores)
                except:
                    logger.info("Can not find \'runConfig.useSysCores, using default Value \'")

                if self.useSysCores:
                    self.MaximumThread=multiprocessing.cpu_count()
                    logger.debug("Auto Detector System CPU_COUNT :{}".format(self.MaximumThread))
                    if self.MaximumThread < 2:
                        self.MaximumThread=2
                        logger.debug("System Cores Does not seems right, using CPU_CORES {}".format(self.MaximumThread))
                else:
                    self.MaximumThread=int((data["runConfig"]["ncore"]))
                    logger.debug("Use User input CPU_COUNT :{}".format(self.MaximumThread))

                self.gapScanMin = data["gapScan"]["min"]
                self.gapScanMax = data["gapScan"]["max"]
                self.gapScanNBin = data["gapScan"]["nbin"]

        else:
            raise IOError("Config file: {} CAN NOT FIND".format(fname))
        # print(data)

        # self.krScanMin=float(data["krScan"]["min"])

    # def HamiSingleScaner(self, SpacialPos=[],gap=1.0):
    #     TaseCal=Tase.Henlarge("topological")
    #     val = la.eigvalsh(TaseCal.HamiltonianMatrix(spacialPos=SpacialPos,gap=gap))
    #
    #     d = val[8] - val[7]
    #     if d < 0.01:
    #         return {"Pos":SpacialPos,"gap":gap,"HamiVal":val.tolist()}
    #     return None


    def HamiSingleScaner(self, SpacialPos=[],gap=1.0):
        TaseCal=Tase.Henlarge("topological")
        val = la.eigvalsh(TaseCal.HamiltonianMatrix(spacialPos=SpacialPos,gap=gap))

        return {"Pos":SpacialPos,"gap":gap,"HamiVal":val.tolist()}
        return None


    def _threadCallBack(self, threadID,totalCount):

        if not 0 in self._progressBarArray.keys():
            self._progressBarArray[0] = tqdm(total=self.MaximumJobs, position=0)
            self._progressBarArray[0].set_description("Total Jobs Process ")
        if self._finishedWorkCount > self._finishedWorkCount_Update:
            self._progressBarArray[0].update(self._finishedWorkCount-self._finishedWorkCount_Update)
            self._finishedWorkCount_Update=self._finishedWorkCount

        if not threadID in self._progressBarArray.keys():
            self._progressBarArray[threadID]=tqdm(total=totalCount,position=threadID+1)
            self._progressBarArray[threadID].set_description("subProcess{:4} ".format(threadID))
        self._progressBarArray[threadID].update(1)
        self._progressBarArray[threadID].refresh()

    def HamiScaner(self, krArray=[],ksArrary=[],kzArrary=[], gapArray =[],_callBack=None):
        ###

        counter=0
        totalCount=len(krArray)*len(ksArrary)
        processID=multiprocessing.Process().name.split('-')[-1]

        ScanResult=[]
        for kr in krArray:
            for ks in ksArrary:
                for kz in kzArrary:
                    for gap in gapArray:
                        SpacialPos=[kr,ks,kz]
                        result = self.HamiSingleScaner(SpacialPos=SpacialPos,gap=gap)
                        if result :
                            ScanResult.append(result)

                self._threadCallBack(int(processID),totalCount)

        self._finishedWorkCount=self._finishedWorkCount+1
        return ScanResult

    def MT_PostProcess(self,result=[]):
        '''
        combine result for multi process thread
        '''
        finalResult={}

        print("Scan report :: {} data point writen to file".format(len(result)))
        for id in range(0,len(result)):
            item=result[id]
            finalResult[id] = item
        self.ResultDumper(result=finalResult)

    def ResultDumper(self,result=[],fname=""):
        '''

        '''
        if not fname:
            fname="./result/dHxy_Result_{}.json".format(datetime.datetime.now().strftime("%Y_%m_%d_%H_%m_%s"))

        with open(fname,"w") as fileio:
            json.dump(result,fileio)
        fileio.close()

    def MTHamiScanner(self):
        '''
        The number of Jobs equals the number of thread that want to process
        Or twice maybe
        '''

        #chop the jobs

        NBinVal=[self.krScanNBin,self.ksScanNBin,self.kzScanNBin, self.gapScanNBin]

        krScanBin = []
        ksScanBin = []
        kzScanBin = []
        gapScanBin= []

        krRange  = np.linspace(self.krScanMin,self.krScanMax,self.krScanNBin)
        ksRange  = np.linspace(self.ksScanMin, self.ksScanMax, self.ksScanNBin)
        kzRange  = np.linspace(self.kzScanMin, self.kzScanMax, self.kzScanNBin)
        gapRange = np.linspace(self.gapScanMin, self.gapScanMax, self.gapScanNBin)

        if self.krScanNBin >= max(NBinVal):
            if max(NBinVal) >= 2*self.MaximumJobs:
                nbin= int(max(NBinVal)/self.MaximumJobs)
                startIndex=0
                endIndex=len(krRange)
                for x in range(0, endIndex,nbin):
                    if x+nbin < endIndex + 1:
                        val=krRange[x:x+nbin]
                        krScanBin.append(val)
                    elif x < endIndex:
                        val = krRange[x:-1]
                        krScanBin.append(val)

            else:
                for item in krRange:
                    krScanBin.append([item])

        elif self.ksScanNBin >= max(NBinVal):
            if max(NBinVal) >= 2*self.MaximumJobs:
                nbin= int(max(NBinVal)/self.MaximumJobs)
                startIndex=0
                endIndex=len(ksRange)
                for x in range(0, endIndex,nbin):
                    if x+nbin < endIndex + 1:
                        val = ksRange[x:x + nbin]
                        ksScanBin.append(val)
                    elif x < endIndex:
                        val = ksRange[x:-1]
                        ksScanBin.append(val)
            else:
                for item in krRange:
                    ksScanBin.append([item])

        elif self.kzScanNBin >= max(NBinVal):
            if max(NBinVal) >= 2 * self.MaximumJobs:
                nbin = int(max(NBinVal) / self.MaximumJobs)
                startIndex = 0
                endIndex = len(kzRange)
                for x in range(0, endIndex, nbin):
                    if x+nbin < endIndex + 1:
                        val = kzRange[x:x + nbin]
                        kzScanBin.append(val)
                    elif x < endIndex:
                        val = kzRange[x:-1]
                        kzScanBin.append(val)
            else:
                for item in krRange:
                    kzScanBin.append([item])
        elif self.gapScanNBin >= max(NBinVal):
            if max(NBinVal) >= 2 * self.MaximumJobs:
                nbin = int(max(NBinVal) / self.MaximumJobs)
                startIndex = 0
                endIndex = len(gapRange)
                for x in range(0, endIndex, nbin):
                    if x+nbin < endIndex + 1:
                        val = gapRange[x:x + nbin]
                        gapScanBin.append(val)
                    elif x < endIndex:
                        val = gapRange[x:-1]
                        gapScanBin.append(val)
            else:
                for item in krRange:
                    gapScanBin.append([item])

        if len(krScanBin) == 0 :
            krScanBin.append(krRange)
        if len(ksScanBin) == 0:
            ksScanBin.append(ksRange)
        if len(kzScanBin) == 0:
            kzScanBin.append(kzRange)
        if len(gapScanBin) == 0:
            gapScanBin.append(gapRange)

        # start thread pool
        logger.debug("kr: {}, ks: {}, kz: {},  gap: {} ".format(len(krScanBin),len(ksScanBin),len(kzScanBin),len(gapScanBin)))

        # multi thread scanner
        pool=ThreadPool(self.MaximumThread)
        threadResult=[]

        for kr in krScanBin:
            for kz in kzScanBin:
                for ks in ksScanBin:
                    for gap in gapScanBin:
                        threadResult.append(pool.apply_async(self.HamiScaner, args=[kr,ks,kz,gap]))

        pool.close()
        pool.join()

        result=[]
        for r in threadResult:
            result=result+r.get()

        self.MT_PostProcess(result=result)

    def tester(self):
        pass


if __name__ == '__main__':
    test=WireFinder()
    if len(sys.argv) == 2:
        if sys.argv[1].endswith('.json'):
            test.LoadRunConfig(fname=sys.argv[1])
        else:
            test.LoadRunConfig()
    else:
        test.LoadRunConfig()
    # test.LoadRunConfig()
    test.MTHamiScanner()


