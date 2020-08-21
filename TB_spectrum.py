import TaSeHamiltonian as TaSe
import numpy as np
from numpy import linalg as la
from math import cos, sin, sqrt, pi
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.linalg import block_diag
from matplotlib.animation import FuncAnimation


#parameters
a, c, v0, u1, u1t, u2 = TaSe.a, TaSe.c, TaSe.v0, TaSe.u1, TaSe.u1t, TaSe.u2

def band2Dkx(HAB,gap,range,ky=0,kz = pi,delta = "dHxy"):
    """
    plot along kx at kz and ky
    :param HAB:Hamiltonian from a class
    :param gap:
    :param range: test.x_range
    :param ky:
    :return:
    """
    numpoints = 200
    kxdata = np.linspace(-range,range, numpoints)
    val = []

    for kx in kxdata:
        val.append(la.eigvalsh(HAB( kx, ky, kz, gap)))

    fig, ax = plt.subplots()

    ax.plot(kxdata, val)
    plt.title(
        str(delta)+": gap=" + str(gap) + "\n" + "  v0=" + str(v0) + "  u1=" + str(u1) + "  u1t=" + str(u1t) + "  u2=" + str(u2) )
    plt.xlabel("k_x")
    plt.show()
    return None

def band2Dkr(HAB,gap,range,ks=0,kz = pi,delta = "dHxy"):
    """
    plot along kr at ks and kz
    :param HAB:Hamiltonian from a class
    :param gap:
    :param range: test.r_range
    :param ky:
    :return:
    """
    numpoints = 200
    krdata = np.linspace(-range,range, numpoints)
    val = []

    for kr in krdata:
        val.append(la.eigvalsh(HAB( kr, ks, kz, gap)))

    fig, ax = plt.subplots()

    ax.plot(krdata, val)
    plt.title(
        str(delta)+": gap=" + str(gap) + "\n" + "  v0=" + str(v0) + "  u1=" + str(u1) + "  u1t=" + str(u1t) + "  u2=" + str(u2) )
    plt.xlabel("k_r")
    plt.show()
    return None





if __name__ == '__main__':
    test = TaSe.Henlarge()
    band2Dkx(test.Hz,0,test.x_range,delta="dHz")
