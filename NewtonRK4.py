import numpy as np
from matplotlib import pyplot as plt

def d2r(r):
    return (-gamma / r ** 2 + l ** 2 / mu / r ** 3) / mu

def phidot(r):
    return (l/mu/r**2)


#Define parameters
G= 6.67e-11*31536000**2 #m^3/kg/year^2
m1 = 5.972e24 #kg
m2 = 1.989e30 #kg
m = m1 + m2
mu = m1*m2/(m1+m2) #kg
r0 = 151.6e9 #m
phidot0 = np.pi*2  #rad/year
l = phidot0*mu*(r0**2)
gamma = G*m1*m2


# Create coordinate lists
r = []
dr = []
phi = []
t = []
r1 = []
r2 = []


# Initial conditions
r.append(r0)
dr.append(0)
phi.append(0)
t.append(0)
r1.append((m2/m*r[0]))
r2.append(m1/m*r[0])


def rk4(t0, tf, h):
    n = (int)((tf-t0)/h)
    i = 0
    while phi[i]<2*np.pi:
        l1 = h * d2r(r[i-1])
        l2 = h * (d2r(r[i-1]+l1/2))
        l3 = h * (d2r(r[i-1]+l2/2))
        l4 = h * (d2r(r[i-1]+l3))
        dr.append(dr[i - 1] + l1/6 + l2/3 +l3/3 +l4/6)

        k1 = h * dr[i - 1]
        k2 = h * (dr[i - 1] + k1/2)
        k3 = h * (dr[i-1] + k2/2)
        k4 = h * (dr[i-1] + k3)
        r.append(r[i-1] + k1/6 + k2/3 +k3/3 +k4/6)

        d2r.append = (-gamma / r[i] ** 2 + l ** 2 / mu / r[i] ** 3) / mu

        j1 = h * phidot(r[i - 1])
        j2 = h * phidot(r[i - 1] + j1/2)
        j3 = h * phidot(r[i-1] + j2/2)
        j4 = h * phidot(r[i-1] + j3)
        phi.append(phi[i-1] + j1/6 + j2/3 + j3/3 +j4/6)

        r1.append(m2 / m * r[i])
        r2.append(-m1 / m * r[i])
        t.append(t[i-1]+h)
        i += 1


rk4(0,1,1e-4)

fig=plt.figure()
ax1 = fig.add_subplot(2,2,1, projection = 'rectilinear',xlabel = "time in years", ylabel ='separation in meters')
ax1.plot(t,r)
ax2 = fig.add_subplot(2,2,2, projection = 'polar', label = 'phi', ylabel = 'distance from center of mass')
ax2.plot(phi,r)
ax3 = fig.add_subplot(2,1,2, projection = 'polar')
ax3.plot(phi, r1)
ax4 = fig.add_subplot(2,1,2, projection = 'polar')
ax4.plot(phi, r2)

plt.show()
