import numpy as np
from matplotlib import pyplot as plt
from eqsofmotion import *

def rk4(t0, tf, h):
    n = int((tf - t0) / h)
    i = 0

    # Create coordinate lists
    r = np.zeros([n])
    rdot = np.zeros([n])
    phi = np.zeros([n])
    phidot = np.zeros([n])
    t = np.zeros([n])
    r1 = np.zeros([n])
    r2 = np.zeros([n])

    # Initial conditions
    r[0] = r0
    rdot[0] = 0
    phi[0] = 0
    phidot[0] = phidot0
    t[0] = t0
    #r1[0] = m2 / m * r[0]
    #r2[0] = -m1 / m * r[0]

    for i in range(0, n-1):
        j1 = rddot(r[i], rdot[i], phidot[i])
        k1 = phiddot(r[i], rdot[i], phidot[i])
        l1 = rdot[i]
        m1 = phidot[i]

        j2 = rddot(r[i] + h * l1 /2, rdot[i] + h * j1 / 2, phidot[i] + h * k1 /2)
        k2 = phiddot(r[i] + h * l1 /2, rdot[i] + h * j1 / 2, phidot[i] + h * k1 /2)
        l2 = rdot[i] + h * j1 /2
        m2 = phidot[i] + h * k1 / 2

        j3 = rddot(r[i] + h * l2 /2, rdot[i] + h * j2 / 2, phidot[i] + h * k2 /2)
        k3 = phiddot(r[i] + h * l2 / 2, rdot[i] + h * j2 / 2, phidot[i] + h * k2 / 2)
        l3 = rdot[i] + h * j2 / 2
        m3 = phidot[i] + h * k2 / 2

        j4 = rddot(r[i] + h * l3, rdot[i] + h * j3, phidot[i] + h * k3)
        k4 = phiddot(r[i] + h * l3, rdot[i] + h * j3, phidot[i] + h * k3)
        l4 = rdot[i] + h * j3
        m4 = phidot[i] + h * k3

        r[i+1] = r[i] + h/6 * (l1 + 2 * l2 + 2 * l3 + l4)
        rdot[i+1] = rdot[i] +h / 6 * (j1 * 2 * j2 +2 * j3 +j4)
        phidot[i+1] = phidot[i] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        phi[i+1] = phi[i] + h / 6 * (m1 + 2 * m2 + 2 * m3 + m4)
        t[i+1] = t[i]+h
        i+=1


    fig = plt.figure()
    ax1 = fig.add_subplot(2, 2, 1, projection='rectilinear', xlabel="time in years", ylabel='separation in meters')
    ax1.plot(t, r)
    ax2 = fig.add_subplot(2, 2, 2, projection='polar', label='phi', ylabel='distance from center of mass')
    ax2.plot(phi, r)
    #ax3= fig.add_subplot(2, 1, 2, projection='polar')
    #ax3.plot(phi, r1)
    #ax3.plot(phi, r2)
    plt.show()



rk4(0, 1, 1e-3)