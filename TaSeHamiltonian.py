import numpy as np
from math import cos, sin, sqrt, pi
from numpy import linalg as la


sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])
sz = np.array([[1,0],[0,-1]])
e = np.identity(2)
sn = sx+sy
snt = sx-sy

#parameters
a, c, v0, u1, u1t, u2 = 1.0, 1.0, 1, 1, 2, 1

mu = np.array([[u1,0],[0,u1t]])
mut = np.array([[u1t,0],[0,u1]])
gamma2 = [e, sx,sy,sz]


sigma = np.array([e, sx, sy, sz])
G16 = [[[[0]* 4 for i in range(4)] for j in range(4)] for k in range(4)]
G8 = [[[0]* 4 for i in range(4)] for j in range(4)]
G4 = [[0]* 4 for i in range(4)]

for i in range(4):
    for j in range(4):
        G4[i][j] =np.kron(sigma[i], sigma[j])

for i in range(4):
    for j in range(4):
        for k in range(4):
            G8[i][j][k] = np.kron(np.kron(sigma[i], sigma[j]), sigma[k])

for i in range(4):
    for j in range(4):
        for k in range(4):
            for l in range(4):
                G16[i][j][k][l] = np.kron(np.kron(np.kron(sigma[i], sigma[j]), sigma[k]), sigma[l])


"""
This program is to describe the band diagram for enlarged unit cell
 with second and third NN hopping with backscattering term in spin z direction 
 basis is sigma tau mu s
"""

class Hamiltonian():
    def __init__(self):
        pass

    def HamiltonianMatrix(self):
        pass


class Henlarge(Hamiltonian):
    """
    Hamiltonian in enlarged unit cell
    top = 1 is the topological gapping along spin z direction
    top = 0 is the trivial gapping in xy plane
    16 by 16 matrix
    """

    hamilType="trivial"
    x_range = pi
    dim = 16
    x_range = pi
    r_range = pi / sqrt(2)
    def __init__(self, type='trivial'):
        '''
        :param type: trivial/topological
        '''

        if type.lower() == 'trivial':
            self.hamilType='trivial'
        else:
            if type.lower() == "topological":
                self.hamilType='topological'


    #KE term


    def HZA1(self, kz):
        return -v0 * (1 + cos(kz)) * G16[0][1][0][0] + v0 * sin(kz) * G16[0][2][0][0]

    def H3(self, kr, ks, kz):
        """
        The H3 KE term
        """
        h = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              (1 / 2) * (u1 + u1t * cos(sqrt(2) * a * kr) +
                         u1t * cos(sqrt(2) * a * ks) + u1 * cos(sqrt(2) * a * (kr + ks)) +
                         u1t * cos(kz) + u1 * cos(sqrt(2) * a * kr) * cos(kz) +
                         u1 * cos(sqrt(2) * a * ks) * cos(kz) +
                         u1t * cos(sqrt(2) * a * (kr + ks)) * cos(kz) +
                         1.0j * u1t * sin(sqrt(2) * a * kr) + 1.0j * u1 * cos(kz) *
                         sin(sqrt(2) * a * kr) + 1.0j * u1t * sin(sqrt(2) * a * ks) +
                         1.0j * u1 * cos(kz) * sin(sqrt(2) * a * ks) +
                         1.0j * u1 * sin(sqrt(2) * a * (kr + ks)) + 1.0j * u1t * cos(kz) *
                         sin(sqrt(2) * a * (kr + ks)) + 1.0j * u1t * sin(kz) +
                         1.0j * u1 * cos(sqrt(2) * a * kr) * sin(kz) + 1.0j * u1 * cos(sqrt(2) * a * ks) *
                         sin(kz) + 1.0j * u1t * cos(sqrt(2) * a * (kr + ks)) * sin(kz) -
                         u1 * sin(sqrt(2) * a * kr) * sin(kz) - u1 * sin(sqrt(2) * a * ks) *
                         sin(kz) - u1t * sin(sqrt(2) * a * (kr + ks)) * sin(kz)), 0, 0,
              0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   (1 / 2) * (u1 + u1t * cos(sqrt(2) * a * kr) +
                              u1t * cos(sqrt(2) * a * ks) + u1 * cos(sqrt(2) * a * (kr + ks)) +
                              u1t * cos(kz) + u1 * cos(sqrt(2) * a * kr) * cos(kz) +
                              u1 * cos(sqrt(2) * a * ks) * cos(kz) +
                              u1t * cos(sqrt(2) * a * (kr + ks)) * cos(kz) +
                              1.0j * u1t * sin(sqrt(2) * a * kr) + 1.0j * u1 * cos(kz) *
                              sin(sqrt(2) * a * kr) + 1.0j * u1t * sin(sqrt(2) * a * ks) +
                              1.0j * u1 * cos(kz) * sin(sqrt(2) * a * ks) +
                              1.0j * u1 * sin(sqrt(2) * a * (kr + ks)) + 1.0j * u1t * cos(kz) *
                              sin(sqrt(2) * a * (kr + ks)) + 1.0j * u1t * sin(kz) +
                              1.0j * u1 * cos(sqrt(2) * a * kr) * sin(kz) + 1.0j * u1 * cos(sqrt(2) * a * ks) *
                              sin(kz) + 1.0j * u1t * cos(sqrt(2) * a * (kr + ks)) * sin(kz) -
                              u1 * sin(sqrt(2) * a * kr) * sin(kz) - u1 * sin(sqrt(2) * a * ks) *
                              sin(kz) - u1t * sin(sqrt(2) * a * (kr + ks)) * sin(kz)), 0,
                   0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        (1 / 2) * (u1 + u1t * cos(sqrt(2) * a * kr) +
                                   u1 * cos(sqrt(2) * a * (kr - ks)) + u1t * cos(sqrt(2) * a * ks) +
                                   u1t * cos(kz) + u1 * cos(sqrt(2) * a * kr) * cos(kz) +
                                   u1t * cos(sqrt(2) * a * (kr - ks)) * cos(kz) +
                                   u1 * cos(sqrt(2) * a * ks) * cos(kz) +
                                   1.0j * u1t * sin(sqrt(2) * a * kr) + 1.0j * u1 * cos(kz) *
                                   sin(sqrt(2) * a * kr) + 1.0j * u1 * sin(sqrt(2) * a * (kr - ks)) +
                                   1.0j * u1t * cos(kz) * sin(sqrt(2) * a * (kr - ks)) -
                                   1.0j * u1t * sin(sqrt(2) * a * ks) - 1.0j * u1 * cos(kz) *
                                   sin(sqrt(2) * a * ks) + 1.0j * u1t * sin(kz) +
                                   1.0j * u1 * cos(sqrt(2) * a * kr) * sin(kz) +
                                   1.0j * u1t * cos(sqrt(2) * a * (kr - ks)) * sin(kz) +
                                   1.0j * u1 * cos(sqrt(2) * a * ks) * sin(kz) - u1 * sin(sqrt(2) * a * kr) *
                                   sin(kz) - u1t * sin(sqrt(2) * a * (kr - ks)) * sin(kz) +
                                   u1 * sin(sqrt(2) * a * ks) * sin(kz)), 0],
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              (1 / 2) * (u1 + u1t * cos(sqrt(2) * a * kr) +
                         u1 * cos(sqrt(2) * a * (kr - ks)) + u1t * cos(sqrt(2) * a * ks) +
                         u1t * cos(kz) + u1 * cos(sqrt(2) * a * kr) * cos(kz) +
                         u1t * cos(sqrt(2) * a * (kr - ks)) * cos(kz) +
                         u1 * cos(sqrt(2) * a * ks) * cos(kz) +
                         1.0j * u1t * sin(sqrt(2) * a * kr) + 1.0j * u1 * cos(kz) *
                         sin(sqrt(2) * a * kr) + 1.0j * u1 * sin(sqrt(2) * a * (kr - ks)) +
                         1.0j * u1t * cos(kz) * sin(sqrt(2) * a * (kr - ks)) -
                         1.0j * u1t * sin(sqrt(2) * a * ks) - 1.0j * u1 * cos(kz) *
                         sin(sqrt(2) * a * ks) + 1.0j * u1t * sin(kz) +
                         1.0j * u1 * cos(sqrt(2) * a * kr) * sin(kz) +
                         1.0j * u1t * cos(sqrt(2) * a * (kr - ks)) * sin(kz) +
                         1.0j * u1 * cos(sqrt(2) * a * ks) * sin(kz) - u1 * sin(sqrt(2) * a * kr) *
                         sin(kz) - u1t * sin(sqrt(2) * a * (kr - ks)) * sin(kz) +
                         u1 * sin(sqrt(2) * a * ks) * sin(kz))], [0, 0, 0, 0, 0, 0, 0,
                                                                  0, (1 / 2) * (u1 + u1t * cos(sqrt(2) * a * kr) +
                                                                                u1t * cos(sqrt(2) * a * ks) + u1 * cos(
                            sqrt(2) * a * (kr + ks)) +
                                                                                u1t * cos(kz) + u1 * cos(
                            sqrt(2) * a * kr) * cos(kz) +
                                                                                u1 * cos(sqrt(2) * a * ks) * cos(kz) +
                                                                                u1t * cos(
                            sqrt(2) * a * (kr + ks)) * cos(kz) +
                                                                                u1 * sin(sqrt(2) * a * kr) * sin(
                            kz) + u1 * sin(sqrt(2) * a * ks) *
                                                                                sin(kz) + u1t * sin(
                            sqrt(2) * a * (kr + ks)) * sin(kz) -
                                                                                1.0j * ((-u1t) * sin(
                            sqrt(2) * a * kr) - u1 * cos(kz) *
                                                                                        sin(sqrt(
                                                                                            2) * a * kr) - u1t * sin(
                                    sqrt(2) * a * ks) -
                                                                                        u1 * cos(kz) * sin(
                                    sqrt(2) * a * ks) -
                                                                                        u1 * sin(
                                    sqrt(2) * a * (kr + ks)) - u1t * cos(kz) *
                                                                                        sin(sqrt(2) * a * (
                                                                                                kr + ks)) + u1t * sin(
                                    kz) +
                                                                                        u1 * cos(
                                    sqrt(2) * a * kr) * sin(kz) + u1 * cos(sqrt(2) * a * ks) *
                                                                                        sin(kz) + u1t * cos(
                                    sqrt(2) * a * (kr + ks)) * sin(kz))), 0,
                                                                  0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                                      (1 / 2) * (u1 + u1t * cos(
                                                                                          sqrt(2) * a * kr) +
                                                                                                 u1t * cos(sqrt(
                                                                                                  2) * a * ks) + u1 * cos(
                                                                                                  sqrt(2) * a * (
                                                                                                          kr + ks)) +
                                                                                                 u1t * cos(
                                                                                                  kz) + u1 * cos(sqrt(
                                                                                                  2) * a * kr) * cos(
                                                                                                  kz) +
                                                                                                 u1 * cos(sqrt(
                                                                                                  2) * a * ks) * cos(
                                                                                                  kz) +
                                                                                                 u1t * cos(
                                                                                                  sqrt(2) * a * (
                                                                                                          kr + ks)) * cos(
                                                                                                  kz) +
                                                                                                 u1 * sin(sqrt(
                                                                                                  2) * a * kr) * sin(
                                                                                                  kz) + u1 * sin(
                                                                                                  sqrt(2) * a * ks) *
                                                                                                 sin(kz) + u1t * sin(
                                                                                                  sqrt(2) * a * (
                                                                                                          kr + ks)) * sin(
                                                                                                  kz) -
                                                                                                 1.0j * ((-u1t) * sin(
                                                                                                  sqrt(
                                                                                                      2) * a * kr) - u1 * cos(
                                                                                                  kz) *
                                                                                                         sin(sqrt(
                                                                                                             2) * a * kr) - u1t * sin(
                                                                                                          sqrt(
                                                                                                              2) * a * ks) -
                                                                                                         u1 * cos(
                                                                                                          kz) * sin(
                                                                                                          sqrt(
                                                                                                              2) * a * ks) -
                                                                                                         u1 * sin(sqrt(
                                                                                                          2) * a * (
                                                                                                                          kr + ks)) - u1t * cos(
                                                                                                          kz) *
                                                                                                         sin(sqrt(
                                                                                                             2) * a * (
                                                                                                                     kr + ks)) + u1t * sin(
                                                                                                          kz) +
                                                                                                         u1 * cos(sqrt(
                                                                                                          2) * a * kr) * sin(
                                                                                                          kz) + u1 * cos(
                                                                                                          sqrt(
                                                                                                              2) * a * ks) *
                                                                                                         sin(
                                                                                                             kz) + u1t * cos(
                                                                                                          sqrt(
                                                                                                              2) * a * (
                                                                                                                  kr + ks)) * sin(
                                                                                                          kz))), 0,
                                                                                      0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              (1 / 2) * (u1 + u1t * cos(sqrt(2) * a * kr) +
                         u1 * cos(sqrt(2) * a * (kr - ks)) + u1t * cos(sqrt(2) * a * ks) +
                         u1t * cos(kz) + u1 * cos(sqrt(2) * a * kr) * cos(kz) +
                         u1t * cos(sqrt(2) * a * (kr - ks)) * cos(kz) +
                         u1 * cos(sqrt(2) * a * ks) * cos(kz) + u1 * sin(sqrt(2) * a * kr) *
                         sin(kz) + u1t * sin(sqrt(2) * a * (kr - ks)) * sin(kz) -
                         u1 * sin(sqrt(2) * a * ks) * sin(kz) -
                         1.0j * ((-u1t) * sin(sqrt(2) * a * kr) - u1 * cos(kz) *
                                 sin(sqrt(2) * a * kr) - u1 * sin(sqrt(2) * a * (kr - ks)) -
                                 u1t * cos(kz) * sin(sqrt(2) * a * (kr - ks)) +
                                 u1t * sin(sqrt(2) * a * ks) + u1 * cos(kz) *
                                 sin(sqrt(2) * a * ks) + u1t * sin(kz) +
                                 u1 * cos(sqrt(2) * a * kr) * sin(kz) +
                                 u1t * cos(sqrt(2) * a * (kr - ks)) * sin(kz) +
                                 u1 * cos(sqrt(2) * a * ks) * sin(kz))), 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              (1 / 2) * (u1 + u1t * cos(sqrt(2) * a * kr) +
                         u1 * cos(sqrt(2) * a * (kr - ks)) + u1t * cos(sqrt(2) * a * ks) +
                         u1t * cos(kz) + u1 * cos(sqrt(2) * a * kr) * cos(kz) +
                         u1t * cos(sqrt(2) * a * (kr - ks)) * cos(kz) +
                         u1 * cos(sqrt(2) * a * ks) * cos(kz) + u1 * sin(sqrt(2) * a * kr) *
                         sin(kz) + u1t * sin(sqrt(2) * a * (kr - ks)) * sin(kz) -
                         u1 * sin(sqrt(2) * a * ks) * sin(kz) -
                         1.0j * ((-u1t) * sin(sqrt(2) * a * kr) - u1 * cos(kz) *
                                 sin(sqrt(2) * a * kr) - u1 * sin(sqrt(2) * a * (kr - ks)) -
                                 u1t * cos(kz) * sin(sqrt(2) * a * (kr - ks)) +
                                 u1t * sin(sqrt(2) * a * ks) + u1 * cos(kz) *
                                 sin(sqrt(2) * a * ks) + u1t * sin(kz) +
                                 u1 * cos(sqrt(2) * a * kr) * sin(kz) +
                                 u1t * cos(sqrt(2) * a * (kr - ks)) * sin(kz) +
                                 u1 * cos(sqrt(2) * a * ks) * sin(kz))), 0, 0, 0, 0],
             [0, 0, 0, 0, (1 / 2) * (u1 + u1t * cos(sqrt(2) * a * kr) +
                                     u1t * cos(sqrt(2) * a * ks) + u1 * cos(sqrt(2) * a * (kr + ks)) +
                                     u1t * cos(kz) + u1 * cos(sqrt(2) * a * kr) * cos(kz) +
                                     u1 * cos(sqrt(2) * a * ks) * cos(kz) +
                                     u1t * cos(sqrt(2) * a * (kr + ks)) * cos(kz) -
                                     1.0j * u1t * sin(sqrt(2) * a * kr) - 1.0j * u1 * cos(kz) *
                                     sin(sqrt(2) * a * kr) - 1.0j * u1t * sin(sqrt(2) * a * ks) -
                                     1.0j * u1 * cos(kz) * sin(sqrt(2) * a * ks) -
                                     1.0j * u1 * sin(sqrt(2) * a * (kr + ks)) - 1.0j * u1t * cos(kz) *
                                     sin(sqrt(2) * a * (kr + ks)) + 1.0j * u1t * sin(kz) +
                                     1.0j * u1 * cos(sqrt(2) * a * kr) * sin(kz) + 1.0j * u1 * cos(sqrt(2) * a * ks) *
                                     sin(kz) + 1.0j * u1t * cos(sqrt(2) * a * (kr + ks)) * sin(kz) +
                                     u1 * sin(sqrt(2) * a * kr) * sin(kz) + u1 * sin(sqrt(2) * a * ks) *
                                     sin(kz) + u1t * sin(sqrt(2) * a * (kr + ks)) * sin(kz)), 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0,
                                           (1 / 2) * (u1 + u1t * cos(sqrt(2) * a * kr) +
                                                      u1t * cos(sqrt(2) * a * ks) + u1 * cos(sqrt(2) * a * (kr + ks)) +
                                                      u1t * cos(kz) + u1 * cos(sqrt(2) * a * kr) * cos(kz) +
                                                      u1 * cos(sqrt(2) * a * ks) * cos(kz) +
                                                      u1t * cos(sqrt(2) * a * (kr + ks)) * cos(kz) -
                                                      1.0j * u1t * sin(sqrt(2) * a * kr) - 1.0j * u1 * cos(kz) *
                                                      sin(sqrt(2) * a * kr) - 1.0j * u1t * sin(sqrt(2) * a * ks) -
                                                      1.0j * u1 * cos(kz) * sin(sqrt(2) * a * ks) -
                                                      1.0j * u1 * sin(sqrt(2) * a * (kr + ks)) - 1.0j * u1t * cos(kz) *
                                                      sin(sqrt(2) * a * (kr + ks)) + 1.0j * u1t * sin(kz) +
                                                      1.0j * u1 * cos(sqrt(2) * a * kr) * sin(kz) + 1.0j * u1 * cos(
                                                       sqrt(2) * a * ks) *
                                                      sin(kz) + 1.0j * u1t * cos(sqrt(2) * a * (kr + ks)) * sin(kz) +
                                                      u1 * sin(sqrt(2) * a * kr) * sin(kz) + u1 * sin(
                                                       sqrt(2) * a * ks) *
                                                      sin(kz) + u1t * sin(sqrt(2) * a * (kr + ks)) * sin(kz)), 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0,
                                                                     (1 / 2) * (u1 + u1t * cos(sqrt(2) * a * kr) +
                                                                                u1 * cos(
                                                                                 sqrt(2) * a * (kr - ks)) + u1t * cos(
                                                                                 sqrt(2) * a * ks) +
                                                                                u1t * cos(kz) + u1 * cos(
                                                                                 sqrt(2) * a * kr) * cos(kz) +
                                                                                u1t * cos(
                                                                                 sqrt(2) * a * (kr - ks)) * cos(kz) +
                                                                                u1 * cos(sqrt(2) * a * ks) * cos(kz) -
                                                                                1.0j * u1t * sin(
                                                                                 sqrt(2) * a * kr) - 1.0j * u1 * cos(
                                                                                 kz) *
                                                                                sin(sqrt(2) * a * kr) - 1.0j * u1 * sin(
                                                                                 sqrt(2) * a * (kr - ks)) -
                                                                                1.0j * u1t * cos(kz) * sin(
                                                                                 sqrt(2) * a * (kr - ks)) +
                                                                                1.0j * u1t * sin(
                                                                                 sqrt(2) * a * ks) + 1.0j * u1 * cos(
                                                                                 kz) *
                                                                                sin(sqrt(
                                                                                    2) * a * ks) + 1.0j * u1t * sin(
                                                                                 kz) +
                                                                                1.0j * u1 * cos(sqrt(2) * a * kr) * sin(
                                                                                 kz) +
                                                                                1.0j * u1t * cos(
                                                                                 sqrt(2) * a * (kr - ks)) * sin(kz) +
                                                                                1.0j * u1 * cos(sqrt(2) * a * ks) * sin(
                                                                                 kz) + u1 * sin(sqrt(2) * a * kr) *
                                                                                sin(kz) + u1t * sin(
                                                                                 sqrt(2) * a * (kr - ks)) * sin(kz) -
                                                                                u1 * sin(sqrt(2) * a * ks) * sin(kz)),
                                                                     0, 0, 0, 0, 0, 0, 0, 0,
                                                                     0], [0, 0, 0, 0, 0, 0, 0,
                                                                          (1 / 2) * (u1 + u1t * cos(sqrt(2) * a * kr) +
                                                                                     u1 * cos(sqrt(2) * a * (
                                                                                          kr - ks)) + u1t * cos(
                                                                                      sqrt(2) * a * ks) +
                                                                                     u1t * cos(kz) + u1 * cos(
                                                                                      sqrt(2) * a * kr) * cos(kz) +
                                                                                     u1t * cos(
                                                                                      sqrt(2) * a * (kr - ks)) * cos(
                                                                                      kz) +
                                                                                     u1 * cos(sqrt(2) * a * ks) * cos(
                                                                                      kz) -
                                                                                     1.0j * u1t * sin(sqrt(
                                                                                      2) * a * kr) - 1.0j * u1 * cos(
                                                                                      kz) *
                                                                                     sin(sqrt(
                                                                                         2) * a * kr) - 1.0j * u1 * sin(
                                                                                      sqrt(2) * a * (kr - ks)) -
                                                                                     1.0j * u1t * cos(kz) * sin(
                                                                                      sqrt(2) * a * (kr - ks)) +
                                                                                     1.0j * u1t * sin(sqrt(
                                                                                      2) * a * ks) + 1.0j * u1 * cos(
                                                                                      kz) *
                                                                                     sin(sqrt(
                                                                                         2) * a * ks) + 1.0j * u1t * sin(
                                                                                      kz) +
                                                                                     1.0j * u1 * cos(
                                                                                      sqrt(2) * a * kr) * sin(kz) +
                                                                                     1.0j * u1t * cos(
                                                                                      sqrt(2) * a * (kr - ks)) * sin(
                                                                                      kz) +
                                                                                     1.0j * u1 * cos(
                                                                                      sqrt(2) * a * ks) * sin(
                                                                                      kz) + u1 * sin(sqrt(2) * a * kr) *
                                                                                     sin(kz) + u1t * sin(
                                                                                      sqrt(2) * a * (kr - ks)) * sin(
                                                                                      kz) -
                                                                                     u1 * sin(sqrt(2) * a * ks) * sin(
                                                                                      kz)), 0, 0, 0, 0, 0, 0, 0, 0],
             [(1 / 2) * (u1 + u1t * cos(sqrt(2) * a * kr) +
                         u1t * cos(sqrt(2) * a * ks) + u1 * cos(sqrt(2) * a * (kr + ks)) +
                         u1t * cos(kz) + u1 * cos(sqrt(2) * a * kr) * cos(kz) +
                         u1 * cos(sqrt(2) * a * ks) * cos(kz) +
                         u1t * cos(sqrt(2) * a * (kr + ks)) * cos(kz) -
                         u1 * sin(sqrt(2) * a * kr) * sin(kz) - u1 * sin(sqrt(2) * a * ks) *
                         sin(kz) - u1t * sin(sqrt(2) * a * (kr + ks)) * sin(kz) -
                         1.0j * (u1t * sin(sqrt(2) * a * kr) + u1 * cos(kz) *
                                 sin(sqrt(2) * a * kr) + u1t * sin(sqrt(2) * a * ks) +
                                 u1 * cos(kz) * sin(sqrt(2) * a * ks) +
                                 u1 * sin(sqrt(2) * a * (kr + ks)) + u1t * cos(kz) *
                                 sin(sqrt(2) * a * (kr + ks)) + u1t * sin(kz) +
                                 u1 * cos(sqrt(2) * a * kr) * sin(kz) + u1 * cos(sqrt(2) * a * ks) *
                                 sin(kz) + u1t * cos(sqrt(2) * a * (kr + ks)) * sin(kz))), 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
             [0, (1 / 2) * (u1 + u1t * cos(sqrt(2) * a * kr) +
                            u1t * cos(sqrt(2) * a * ks) + u1 * cos(sqrt(2) * a * (kr + ks)) +
                            u1t * cos(kz) + u1 * cos(sqrt(2) * a * kr) * cos(kz) +
                            u1 * cos(sqrt(2) * a * ks) * cos(kz) +
                            u1t * cos(sqrt(2) * a * (kr + ks)) * cos(kz) -
                            u1 * sin(sqrt(2) * a * kr) * sin(kz) - u1 * sin(sqrt(2) * a * ks) *
                            sin(kz) - u1t * sin(sqrt(2) * a * (kr + ks)) * sin(kz) -
                            1.0j * (u1t * sin(sqrt(2) * a * kr) + u1 * cos(kz) *
                                    sin(sqrt(2) * a * kr) + u1t * sin(sqrt(2) * a * ks) +
                                    u1 * cos(kz) * sin(sqrt(2) * a * ks) +
                                    u1 * sin(sqrt(2) * a * (kr + ks)) + u1t * cos(kz) *
                                    sin(sqrt(2) * a * (kr + ks)) + u1t * sin(kz) +
                                    u1 * cos(sqrt(2) * a * kr) * sin(kz) + u1 * cos(sqrt(2) * a * ks) *
                                    sin(kz) + u1t * cos(sqrt(2) * a * (kr + ks)) * sin(kz))), 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, (1 / 2) * (u1 + u1t * cos(sqrt(2) * a * kr) +
                               u1 * cos(sqrt(2) * a * (kr - ks)) + u1t * cos(sqrt(2) * a * ks) +
                               u1t * cos(kz) + u1 * cos(sqrt(2) * a * kr) * cos(kz) +
                               u1t * cos(sqrt(2) * a * (kr - ks)) * cos(kz) +
                               u1 * cos(sqrt(2) * a * ks) * cos(kz) - u1 * sin(sqrt(2) * a * kr) *
                               sin(kz) - u1t * sin(sqrt(2) * a * (kr - ks)) * sin(kz) +
                               u1 * sin(sqrt(2) * a * ks) * sin(kz) -
                               1.0j * (u1t * sin(sqrt(2) * a * kr) + u1 * cos(kz) *
                                       sin(sqrt(2) * a * kr) + u1 * sin(sqrt(2) * a * (kr - ks)) +
                                       u1t * cos(kz) * sin(sqrt(2) * a * (kr - ks)) -
                                       u1t * sin(sqrt(2) * a * ks) - u1 * cos(kz) *
                                       sin(sqrt(2) * a * ks) + u1t * sin(kz) +
                                       u1 * cos(sqrt(2) * a * kr) * sin(kz) +
                                       u1t * cos(sqrt(2) * a * (kr - ks)) * sin(kz) +
                                       u1 * cos(sqrt(2) * a * ks) * sin(kz))), 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0], [0, 0, 0,
                                  (1 / 2) * (u1 + u1t * cos(sqrt(2) * a * kr) +
                                             u1 * cos(sqrt(2) * a * (kr - ks)) + u1t * cos(sqrt(2) * a * ks) +
                                             u1t * cos(kz) + u1 * cos(sqrt(2) * a * kr) * cos(kz) +
                                             u1t * cos(sqrt(2) * a * (kr - ks)) * cos(kz) +
                                             u1 * cos(sqrt(2) * a * ks) * cos(kz) - u1 * sin(sqrt(2) * a * kr) *
                                             sin(kz) - u1t * sin(sqrt(2) * a * (kr - ks)) * sin(kz) +
                                             u1 * sin(sqrt(2) * a * ks) * sin(kz) -
                                             1.0j * (u1t * sin(sqrt(2) * a * kr) + u1 * cos(kz) *
                                                     sin(sqrt(2) * a * kr) + u1 * sin(sqrt(2) * a * (kr - ks)) +
                                                     u1t * cos(kz) * sin(sqrt(2) * a * (kr - ks)) -
                                                     u1t * sin(sqrt(2) * a * ks) - u1 * cos(kz) *
                                                     sin(sqrt(2) * a * ks) + u1t * sin(kz) +
                                                     u1 * cos(sqrt(2) * a * kr) * sin(kz) +
                                                     u1t * cos(sqrt(2) * a * (kr - ks)) * sin(kz) +
                                                     u1 * cos(sqrt(2) * a * ks) * sin(kz))), 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0]]
        h = np.array(h)
        return h

    def HZA4(self, kr, ks):
        return 0.5 * u2 * (cos(ks * sqrt(2) * a) - cos(kr * sqrt(2) * a)) * G16[0][3][3][0]

    def dHtopz(self, kr, ks):
        return 0.5*np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0,
  1 + cos((a* kr)/sqrt(2) + (a* ks)/sqrt(2)) -
   1.0j*sin((a* kr)/sqrt(2) + (a* ks)/sqrt(2)), 0, 0, 0,
  1.0j - 1.0j*cos((a* kr)/sqrt(2) + (a* ks)/sqrt(2)) -
   sin((a* kr)/sqrt(2) + (a* ks)/sqrt(2)), 0, 0], [0, 0, 0, 0, 0, 0, 0,
  0, -1 - cos((a* kr)/sqrt(2) + (a* ks)/sqrt(2)) +
   1.0j*sin((a* kr)/sqrt(2) + (a* ks)/sqrt(2)), 0, 0, 0,
  1.0j - 1.0j*cos((a* kr)/sqrt(2) + (a* ks)/sqrt(2)) -
   sin((a* kr)/sqrt(2) + (a* ks)/sqrt(2)), 0, 0, 0], [0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0,
  1.0j + 1.0j*cos((a* kr)/sqrt(2) - (a* ks)/sqrt(2)) +
   sin((a* kr)/sqrt(2) - (a* ks)/sqrt(2)), 0, 0,
  0, -1 + cos((a* kr)/sqrt(2) - (a* ks)/sqrt(2)) -
   1.0j*sin((a* kr)/sqrt(2) - (a* ks)/sqrt(2))], [0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 1.0j + 1.0j*cos((a* kr)/sqrt(2) - (a* ks)/sqrt(2)) +
   sin((a* kr)/sqrt(2) - (a* ks)/sqrt(2)), 0, 0, 0,
  1 - cos((a* kr)/sqrt(2) - (a* ks)/sqrt(2)) +
   1.0j*sin((a* kr)/sqrt(2) - (a* ks)/sqrt(2)), 0], [0, 0, 0, 0, 0, 0, 0,
  0, 0, 1.0j - 1.0j*cos((a* kr)/sqrt(2) + (a* ks)/sqrt(2)) -
   sin((a* kr)/sqrt(2) + (a* ks)/sqrt(2)), 0, 0,
  0, -1 - cos((a* kr)/sqrt(2) + (a* ks)/sqrt(2)) +
   1.0j*sin((a* kr)/sqrt(2) + (a* ks)/sqrt(2)), 0, 0], [0, 0, 0, 0, 0, 0,
  0, 0, 1.0j - 1.0j*cos((a* kr)/sqrt(2) + (a* ks)/sqrt(2)) -
   sin((a* kr)/sqrt(2) + (a* ks)/sqrt(2)), 0, 0, 0,
  1 + cos((a* kr)/sqrt(2) + (a* ks)/sqrt(2)) -
   1.0j*sin((a* kr)/sqrt(2) + (a* ks)/sqrt(2)), 0, 0, 0], [0, 0, 0, 0, 0,
  0, 0, 0, 0, 0,
  0, -1 + cos((a* kr)/sqrt(2) - (a* ks)/sqrt(2)) -
   1.0j*sin((a* kr)/sqrt(2) - (a* ks)/sqrt(2)), 0, 0,
  0, -1.0j - 1.0j*cos((a* kr)/sqrt(2) - (a* ks)/sqrt(2)) -
   sin((a* kr)/sqrt(2) - (a* ks)/sqrt(2))], [0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 1 - cos((a* kr)/sqrt(2) - (a* ks)/sqrt(2)) +
   1.0j*sin((a* kr)/sqrt(2) - (a* ks)/sqrt(2)), 0, 0,
  0, -1.0j - 1.0j*cos((a* kr)/sqrt(2) - (a* ks)/sqrt(2)) -
   sin((a* kr)/sqrt(2) - (a* ks)/sqrt(2)),
  0], [0, -1 - cos((a* (kr + ks))/sqrt(2)) -
   1.0j*sin((a* (kr + ks))/sqrt(2)), 0, 0,
  0, -1.0j + 1.0j*cos((a* (kr + ks))/sqrt(2)) - sin((a* (kr + ks))/sqrt(2)),
  0, 0, 0, 0, 0, 0, 0, 0, 0,
  0], [1 + cos((a* (kr + ks))/sqrt(2)) + 1.0j*sin((a* (kr + ks))/sqrt(2)),
  0, 0, 0, -1.0j + 1.0j*cos((a* (kr + ks))/sqrt(2)) -
   sin((a* (kr + ks))/sqrt(2)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0,
  0, 0, -1.0j - 1.0j*cos((a* kr)/sqrt(2) - (a* ks)/sqrt(2)) +
   sin((a* kr)/sqrt(2) - (a* ks)/sqrt(2)), 0, 0, 0,
  1 - cos((a* kr)/sqrt(2) - (a* ks)/sqrt(2)) -
   1.0j*sin((a* kr)/sqrt(2) - (a* ks)/sqrt(2)), 0, 0, 0, 0, 0, 0, 0,
  0], [0, 0, -1.0j - 1.0j*cos((a* kr)/sqrt(2) - (a* ks)/sqrt(2)) +
   sin((a* kr)/sqrt(2) - (a* ks)/sqrt(2)), 0, 0,
  0, -1 + cos((a* kr)/sqrt(2) - (a* ks)/sqrt(2)) +
   1.0j*sin((a* kr)/sqrt(2) - (a* ks)/sqrt(2)), 0, 0, 0, 0, 0, 0, 0, 0,
  0], [0, -1.0j + 1.0j*cos((a* (kr + ks))/sqrt(2)) -
   sin((a* (kr + ks))/sqrt(2)), 0, 0, 0,
  1 + cos((a* (kr + ks))/sqrt(2)) + 1.0j*sin((a* (kr + ks))/sqrt(2)), 0, 0,
   0, 0, 0, 0, 0, 0, 0,
  0], [-1.0j + 1.0j*cos((a* (kr + ks))/sqrt(2)) - sin((a* (kr + ks))/sqrt(2)),
   0, 0, 0, -1 - cos((a* (kr + ks))/sqrt(2)) -
   1.0j*sin((a* (kr + ks))/sqrt(2)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0,
   0, 0, 1 - cos((a* kr)/sqrt(2) - (a* ks)/sqrt(2)) -
   1.0j*sin((a* kr)/sqrt(2) - (a* ks)/sqrt(2)), 0, 0, 0,
  1.0j + 1.0j*cos((a* kr)/sqrt(2) - (a* ks)/sqrt(2)) -
   sin((a* kr)/sqrt(2) - (a* ks)/sqrt(2)), 0, 0, 0, 0, 0, 0, 0, 0], [0,
  0, -1 + cos((a* kr)/sqrt(2) - (a* ks)/sqrt(2)) +
   1.0j*sin((a* kr)/sqrt(2) - (a* ks)/sqrt(2)), 0, 0, 0,
  1.0j + 1.0j*cos((a* kr)/sqrt(2) - (a* ks)/sqrt(2)) -
   sin((a* kr)/sqrt(2) - (a* ks)/sqrt(2)), 0, 0, 0, 0, 0, 0, 0, 0, 0]])

    def dHtri(self):
        """
        trivial gapping term
        """
        return G16[1][3][0][0]

    def Hz(self,kr, ks, kz, gap):
        """
        The topological Hamiltonian with spin z gapping
        The total Hamiltonian in kr,ks basis. The intra and inter mass terms are included
        """
        return self.HZA1(kz) + self.H3(kr, ks, kz) + self.HZA4(kr, ks) + gap * self.dHtopz(kr, ks)

    def Htri(self,kr, ks, kz, gap):
        """
        The trivial Hamiltonian in kr,ks coordinates
        :param kr:
        :param ks:
        :param kz:
        :param gap:
        :return: Htrivial
        """
        return self.HZA1(kz) + self.H3(kr, ks, kz) + self.HZA4(kr, ks) + gap  * self.dHtri()

    def HamiltonianMatrix(self,spacialPos=[0,0,0],gap=0):
        if len(spacialPos) != 3:
            raise TypeError

        if self.hamilType == "trivial":
            return self.Htri(kr=spacialPos[0], ks=spacialPos[1], kz=spacialPos[2], gap=gap)
        elif self.hamilType == "topological":
            return self.Hz(kr=spacialPos[0], ks=spacialPos[1], kz=spacialPos[2], gap=gap)


class Horiginal(Hamiltonian):
    '''
    Hamiltonian in original unit cell
    8 by 8 matrix
    '''

    dim = 8
    x_range = pi
    r_range = sqrt(2) * pi

    def __init__(self):
        pass

    # KE term
    def HZA1(self, kz):
        return -v0 * (1 + cos(kz)) * G8[1][0][0] + v0 * sin(kz) * G8[2][0][0]

    def HZA4(self, kx, ky):
        return 0.5 * u2 * (cos((kx + ky) * a) - cos((kx - ky) * a)) * G8[3][3][0]

    def H3(self, kx, ky, kz):
        h = np.kron((cos(kx * a) * sx + cos(ky * a) * cos(kz) * sx - cos(ky * a) * sin(kz) * sy), mu)
        h += np.kron((cos(ky * a) * sx + cos(kx * a) * cos(kz) * sx - cos(kx * a) * sin(kz) * sy), mut)
        return np.kron(h, e)

    def dHtopxy(self, kx, ky):
        g = (cos(kx * a) + cos(ky * a)) * G8[3][0][0] + (cos(kx * a) - cos(ky * a)) * G8[3][3][0]
        g += sin(kx * a) * G8[1][0][1] + sin(ky * a) * G8[1][0][2]
        g += sin(kx * a) * G8[1][3][1] - sin(ky * a) * G8[1][3][2]
        return 0.5 * g

    def Hxy(self, kr, ks, kz, gap):
        """
        The topological Hamiltonian with spin in xy plane.
        The intra and inter mass terms are included
        """
        kx = 1.0/sqrt(2)*(kr+ks)
        ky = 1.0/sqrt(2)*(-kr+ks)
        return self.HZA1(kz) + self.H3(kx, ky, kz) + self.HZA4(kx, ky) + gap * self.dHtopxy(kx, ky)

    def HamiltonianMatrix(self, spacialPos=[0, 0, 0], gap=0):
        if len(spacialPos) != 3:
            raise TypeError
        return self.Hxy(kr=spacialPos[0], ks=spacialPos[1], kz=spacialPos[2], gap=gap)

if __name__ == '__main__':
    print(la.eigvalsh(Henlarge(type='topological').HamiltonianMatrix()))