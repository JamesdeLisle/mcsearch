from math import cos,sin,exp, pow
import cmath as cm
from numpy import real, dot, absolute, array, conjugate, linalg

def getHam(coefficients,p):
    
   
    p0 = p
    mu0	    = coefficients['mu0'].value
    mu1	    = coefficients['mu1'].value
    t0	    = coefficients['t0'].value
    Delta0  = coefficients['Delta0'].value
    t1	    = coefficients['t1'].value
    Delta1  = coefficients['Delta1'].value
    t2	    = coefficients['t2'].value
    Delta2  = coefficients['Delta2'].value
    t3	    = coefficients['t3'].value
    Delta3  = coefficients['Delta3'].value
    t4	    = coefficients['t4'].value
    Delta4  = coefficients['Delta4'].value
      
 
    I = 1j
    
    Ham = 1j*array([[mu0 + t1*cm.exp(I*p0)/2 + cm.exp(-I*p0)*conjugate(t1)/2, t0/2 + t2*cm.exp(I*p0)/2 + cm.exp(-I*p0)*conjugate(t3)/2, -Delta1*cm.exp(I*p0)/2 + cm.exp(-I*p0)*conjugate(Delta1)/2, -Delta0/2 - Delta2*cm.exp(I*p0)/2 + cm.exp(-I*p0)*conjugate(Delta3)/2], [t3*cm.exp(I*p0)/2 + conjugate(t0)/2 + cm.exp(-I*p0)*conjugate(t2)/2, mu1 + t4*cm.exp(I*p0)/2 + cm.exp(-I*p0)*conjugate(t4)/2, -Delta3*cm.exp(I*p0)/2 + conjugate(Delta0)/2 + cm.exp(-I*p0)*conjugate(Delta2)/2, -Delta4*cm.exp(I*p0)/2 + cm.exp(-I*p0)*conjugate(Delta4)/2], [-Delta1*cm.exp(-I*p0)/2 + cm.exp(I*p0)*conjugate(Delta1)/2, -Delta3*cm.exp(-I*p0)/2 + cm.exp(I*p0)*conjugate(Delta2)/2 + conjugate(Delta0)/2, -mu0 - t1*cm.exp(I*p0)/2 - cm.exp(-I*p0)*conjugate(t1)/2, -t0/2 - t2*cm.exp(I*p0)/2 - cm.exp(-I*p0)*conjugate(t3)/2], [-Delta0/2 - Delta2*cm.exp(-I*p0)/2 + cm.exp(I*p0)*conjugate(Delta3)/2, -Delta4*cm.exp(-I*p0)/2 + cm.exp(I*p0)*conjugate(Delta4)/2, -t3*cm.exp(I*p0)/2 - conjugate(t0)/2 - cm.exp(-I*p0)*conjugate(t2)/2, -mu1 - t4*cm.exp(I*p0)/2 - cm.exp(-I*p0)*conjugate(t4)/2]])

    #Ham = array([[mu0 + t1*cos(p0), t0/2.0 + t2*cm.exp(I*p0)/2.0 + t3*cm.exp(-I*p0)/2.0,-I*Delta1*sin(p0), -Delta0/2.0 - Delta2*cm.exp(I*p0)/2 + Delta3*cm.exp(-I*p0)/2.0],[t0/2.0 + t2*cm.exp(-I*p0)/2.0 + t3*cm.exp(I*p0)/2.0, mu1 + t4*cos(p0), Delta0/2.0 + Delta2*cm.exp(-I*p0)/2.0 - Delta3*cm.exp(I*p0)/2.0, -I*Delta4*sin(p0)], [I*Delta1*sin(p0),Delta0/2.0 + Delta2*cm.exp(I*p0)/2.0 - Delta3*cm.exp(-I*p0)/2.0, -mu0 - t1*cos(p0),-t0/2.0 - t2*cm.exp(I*p0)/2 - t3*cm.exp(-I*p0)/2.0], [-Delta0/2.0 - Delta2*cm.exp(-I*p0)/2.0 + Delta3*cm.exp(I*p0)/2.0, I*Delta4*sin(p0), -t0/2.0 - t2*cm.exp(-I*p0)/2.0 - t3*cm.exp(I*p0)/2.0,-mu1 - t4*cos(p0)]])
    
    return Ham

