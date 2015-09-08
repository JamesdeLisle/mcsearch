from math import cos,sin,exp, pow
import cmath as cm
from numpy import real, dot, absolute, array, conjugate, linalg

def getHam(Lambda,p):

    p0 = p
    mu0	    = Lambda[0][0] 
    mu1	    = Lambda[0][1] 
    t0	    = 1j*Lambda[1][2]
    Delta0  = 1j*Lambda[1][3]
    t1	    = Lambda[0][4] 
    Delta1  = Lambda[0][5] 
    t2	    = 1j*Lambda[1][6]
    Delta2  = 1j*Lambda[1][7]
    t3	    = 1j*Lambda[1][8]
    Delta3  = 1j*Lambda[1][9]
    t4	    = Lambda[0][10] 
    Delta4  = Lambda[0][11]  
    I = 1j
    
    Ham = 1j*array([[mu0 + t1*cm.exp(I*p0)/2 + cm.exp(-I*p0)*conjugate(t1)/2, t0/2 + t2*cm.exp(I*p0)/2 + cm.exp(-I*p0)*conjugate(t3)/2, -Delta1*cm.exp(I*p0)/2 + cm.exp(-I*p0)*conjugate(Delta1)/2, -Delta0/2 - Delta2*cm.exp(I*p0)/2 + cm.exp(-I*p0)*conjugate(Delta3)/2], [t3*cm.exp(I*p0)/2 + conjugate(t0)/2 + cm.exp(-I*p0)*conjugate(t2)/2, mu1 + t4*cm.exp(I*p0)/2 + cm.exp(-I*p0)*conjugate(t4)/2, -Delta3*cm.exp(I*p0)/2 + conjugate(Delta0)/2 + cm.exp(-I*p0)*conjugate(Delta2)/2, -Delta4*cm.exp(I*p0)/2 + cm.exp(-I*p0)*conjugate(Delta4)/2], [-Delta1*cm.exp(-I*p0)/2 + cm.exp(I*p0)*conjugate(Delta1)/2, -Delta3*cm.exp(-I*p0)/2 + cm.exp(I*p0)*conjugate(Delta2)/2 + conjugate(Delta0)/2, -mu0 - t1*cm.exp(I*p0)/2 - cm.exp(-I*p0)*conjugate(t1)/2, -t0/2 - t2*cm.exp(I*p0)/2 - cm.exp(-I*p0)*conjugate(t3)/2], [-Delta0/2 - Delta2*cm.exp(-I*p0)/2 + cm.exp(I*p0)*conjugate(Delta3)/2, -Delta4*cm.exp(-I*p0)/2 + cm.exp(I*p0)*conjugate(Delta4)/2, -t3*cm.exp(I*p0)/2 - conjugate(t0)/2 - cm.exp(-I*p0)*conjugate(t2)/2, -mu1 - t4*cm.exp(I*p0)/2 - cm.exp(-I*p0)*conjugate(t4)/2]])

    #Ham = array([[mu0 + t1*cos(p0), t0/2.0 + t2*cm.exp(I*p0)/2.0 + t3*cm.exp(-I*p0)/2.0,-I*Delta1*sin(p0), -Delta0/2.0 - Delta2*cm.exp(I*p0)/2 + Delta3*cm.exp(-I*p0)/2.0],[t0/2.0 + t2*cm.exp(-I*p0)/2.0 + t3*cm.exp(I*p0)/2.0, mu1 + t4*cos(p0), Delta0/2.0 + Delta2*cm.exp(-I*p0)/2.0 - Delta3*cm.exp(I*p0)/2.0, -I*Delta4*sin(p0)], [I*Delta1*sin(p0),Delta0/2.0 + Delta2*cm.exp(I*p0)/2.0 - Delta3*cm.exp(-I*p0)/2.0, -mu0 - t1*cos(p0),-t0/2.0 - t2*cm.exp(I*p0)/2 - t3*cm.exp(-I*p0)/2.0], [-Delta0/2.0 - Delta2*cm.exp(-I*p0)/2.0 + Delta3*cm.exp(I*p0)/2.0, I*Delta4*sin(p0), -t0/2.0 - t2*cm.exp(-I*p0)/2.0 - t3*cm.exp(I*p0)/2.0,-mu1 - t4*cos(p0)]])
    
    return Ham

