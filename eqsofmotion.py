import numpy as np

#Define parameters
G= 6.67e-11*(31536000**2) #m^3/kg/year^2
c = 299792458 * 31536000 #m/yearmu = m1*m2/m #kg
m1 = 5.972e24 #kg
m2 = 1.989e30 #kg
m = m1 + m2
r0 = 151.6e9 #m
phidot0 = np.pi*2  #rad/year
gamma = G*m1*m2
mu = m1*m2/(m)
nu = mu/m
l = phidot0*mu*(r0**2)

def phidot(r):
    return (l/mu/r**2)

def rddot(r, rdot, phidot):
    return -G*m/(r**2) + r * (phidot**2)

def phiddot(r, rdot,phidot):
    return -2*rdot*phidot/r

def a2_5(r, rdot,phidot):
    return (1 * (-24/5*nu*rdot*m/r*(rdot**2 +r**2*phidot**2) -136/15*nu**2*m**2*rdot/r**2))

def b2_5(r,rdot, phidot):
    return 1 * (8*nu*m/5/r*(rdot**2+r**2*phidot**2) + 24/5*nu**2*m**2/r**2)

def a(r, rdot, phidot):
    #return a2_5(r, rdot, phidot)
    return 0

def b(r, rdot, phidot):
    #return b2_5(r, rdot, phidot)
    return 0
