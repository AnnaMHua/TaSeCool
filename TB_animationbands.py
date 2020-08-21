import numpy as np
from numpy import linalg as la
from math import sqrt, pi
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.linalg import block_diag
from matplotlib.animation import FuncAnimation
import TaSeHamiltonian as TaSe


class bandAnimation2D(object):
    '''

    '''

    fig = plt.figure()
    ax = plt.axes(xlim=(-pi, pi), ylim=(-6, 6))
    line, = ax.plot([], [], lw=1)
    lines = [plt.plot([], [])[0] for _ in range(16)]
    plt.title( "kz=0 ,"+"dHz: v0="+str(v0)+", u1="+str(u1)+", u1t="+str(u1t)+", u2="+str(u2))
    wrt= ax.text(0.5,0.8,"gap=0")

    def __init__(self):
        self.ax.set_xlabel("kx")
        pass

    def Get32D(self,a=TaSe.Henlarge,arg=[]):
        a=TaSe.Henlarge("tr")
        a.HamiltonMatrix()
        result=a.Hum(arg[1],arg[0])

    def animatekx(self,i):
        numpoints = 200
        kxdata = np.linspace(-pi,pi,numpoints)
        val = []
        self.wrt.set_text("gap="+str(i*0.1))
        for kx in kxdata:
            val.append(la.eigvalsh(HAB(kx,0, 0, 0.1*i)))
        val = np.array(val)
        for ith in range(8):
            eigen=[item[ith] for item in val]
            self.lines[ith].set_data(kxdata,eigen)
        return tuple(self.lines)+(self.wrt,)

    def animatekr(self,i):
        numpoints = 200
        krdata = np.linspace(-pi/sqrt(2),pi/sqrt(2),numpoints)
        val = []
        self.wrt.set_text("gap=" + str(i * 0.1))
        for kr in krdata:
            val.append(la.eigvalsh(HAB(kr,0, 0,0.1*i)))
        val = np.array(val)
        for ith in range(8):
            eigen=[item[ith] for item in val]
            self.lines[ith].set_data(krdata,eigen)
        return tuple(self.lines)+(self.wrt,)

    def animatekz(self,i):

        numpoints = 200
        kzdata = np.linspace(-pi,pi,numpoints)
        val = []
        self.wrt.set_text("gap=" + str(i * 0.1))
        for kz in kzdata:
            val.append(la.eigvalsh(HAB(0,0, kz,0.1*i)))
        val = np.array(val)
        for ith in range(8):
            eigen=[item[ith] for item in val]
            self.lines[ith].set_data(kzdata,eigen)
        return tuple(self.lines)+(self.wrt,)


    def _dataGenerator(self):
        pass

    def GetAnimation(self, animateTerm="kx", saveName=''):
        '''
        animateTerm : kx, kr, kz
        '''
        animateTerm.lower()
        if 'kx' in animateTerm:
            anim = FuncAnimation(self.fig, self.animatekx, frames=np.linspace(0,30,60), interval=20,blit=True)

        if 'kr' in animateTerm:
            anim = FuncAnimation(self.fig, self.animatekr, frames=np.linspace(0, 70, 70), interval=20, blit=True)

        if 'kz' in animateTerm:
            anim = FuncAnimation(self.fig, self.animatekz, frames=np.linspace(0, 3, 60), interval=20, blit=True)
        plt.show()
        if saveName:
            anim.save(saveName)


if __name__ == '__main__':
    H=TaSe.Henlarge()
    a=bandAnimation2D()
    a.cal

