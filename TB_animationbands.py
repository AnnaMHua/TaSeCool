import numpy as np
from numpy import linalg as la
from math import sqrt, pi
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.linalg import block_diag
from matplotlib.animation import FuncAnimation
import TaSeHamiltonian as TaSe


#parameters
a, c, v0, u1, u1t, u2 = TaSe.a, TaSe.c, TaSe.v0, TaSe.u1, TaSe.u1t, TaSe.u2

class bandAnimation2D(object):
    '''

    '''

    fig = plt.figure()
    ax = plt.axes(xlim=(-pi, pi), ylim=(-40, 40))
    line, = ax.plot([], [], lw=1)
    lines = [plt.plot([], [])[0] for _ in range(16)]

    wrt= ax.text(0.5,0.8,"gap=0")

    def __init__(self,HAB=TaSe.Horiginal()):
        self.ax.set_xlabel("kx")
        self.HAB = HAB
        self.dim = HAB.dim
        self.xrange = HAB.x_range
        self.rrange = HAB.r_range

        self.kz = pi
        self.ks = 0
        self.ky = 0
        pass

    def animatekx(self,i):
        numpoints = 200
        kxdata = np.linspace(-self.xrange,self.xrange,numpoints)
        val = []
        self.wrt.set_text("gap="+str(i))
        for kx in kxdata:
            val.append(la.eigvalsh(self.HAB.HamiltonianMatrix([1.0/sqrt(2)*(kx - self.ky), 1.0/sqrt(2)*(kx + self.ky), self.kz], i)))
        val = np.array(val)
        for ith in range(self.dim):
            eigen=[item[ith] for item in val]
            self.lines[ith].set_data(kxdata,eigen)
        return tuple(self.lines)+(self.wrt,)

    def animatekr(self,i):
        numpoints = 200
        krdata = np.linspace(-self.rrange,self.rrange,numpoints)
        val = []
        self.wrt.set_text("gap=" + str(i ))
        for kr in krdata:
            val.append(la.eigvalsh(self.HAB.HamiltonianMatrix([kr,0, self.kz] ,i)))
        val = np.array(val)
        for ith in range(self.dim):
            eigen=[item[ith] for item in val]
            self.lines[ith].set_data(krdata,eigen)
        return tuple(self.lines)+(self.wrt,)

    def animatekz(self,i):

        numpoints = 200
        kzdata = np.linspace(-pi,pi,numpoints)
        val = []
        self.wrt.set_text("gap=" + str(i))
        for kz in kzdata:
            val.append(la.eigvalsh(self.HAB.HamiltonianMatrix([0,0, self.kz],i)))
        val = np.array(val)
        for ith in range(self.dim):
            eigen=[item[ith] for item in val]
            self.lines[ith].set_data(kzdata,eigen)
        return tuple(self.lines)+(self.wrt,)



    def GetAnimation(self, animateTerm="kx", saveName=''):
        '''
        animateTerm : kx, kr, kz
        '''
        animateTerm.lower()
        if 'kx' in animateTerm:
            anim = FuncAnimation(self.fig, self.animatekx, frames=np.linspace(0,30,60), interval=20,blit=True)

        if 'kr' in animateTerm:
            anim = FuncAnimation(self.fig, self.animatekr, frames=np.linspace(0, 20, 70), interval=20, blit=True)

        if 'kz' in animateTerm:
            anim = FuncAnimation(self.fig, self.animatekz, frames=np.linspace(0, 30, 60), interval=20, blit=True)
        plt.title("kz="+str(self.kz) + ",dHz: v0=" + str(v0) + ", u1=" + str(u1) + ", u1t=" + str(u1t) + ", u2=" + str(u2))
        plt.show()

        if saveName:
            anim.save(saveName)


if __name__ == '__main__':
    test = TaSe.Henlarge("topological")
    print(test.hamilType)
    animate = bandAnimation2D(test)
    animate.GetAnimation(animateTerm="kx")

