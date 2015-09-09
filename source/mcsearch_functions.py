from numpy import real, dot, absolute, array, conjugate, linalg, set_printoptions, sqrt, trace, divide
from generate_hamiltonian import *
import random

def linspace(pMin,pMax,pInt):
    
    space = []
    space.append(pMin)
    inc = (pMax-pMin)/pInt
    for tik in range(0,pInt):
        space.append(space[tik]+inc)
    
    return space

def getCostFunction(coefficients,TRMatrix,momentum_discretisation):
     
    pi = 3.14159
    momentumSpace = linspace(0,pi,momentum_discretisation)
    costFunction = 0.0      
    
    for tik1 in range(0,len(momentumSpace)):
        Ham_PlusP = getHam(coefficients,momentumSpace[tik1])
        Ham_MinusP = getHam(coefficients,-momentumSpace[tik1])
        A = TRMatrix.dot(conjugate(Ham_PlusP).dot(conjugate(TRMatrix).transpose())) - Ham_MinusP     
        if (costFunction < sqrt(trace(A.dot(conjugate(A).transpose()))).real): 
            costFunction = sqrt(trace(A.dot(conjugate(A).transpose()))).real    
    
    return costFunction

def updateStandardDeviation(freqAccepted,standardDeviation):
    
    standardDeviation_temp = standardDeviation
    freqAccepted_ideal = 0.5 
    scalingFactor = absolute(freqAccepted-freqAccepted_ideal)
    
    if freqAccepted > freqAccepted_ideal:
        standardDeviation_temp = standardDeviation/scalingFactor
    else:
        standardDeviation_temp = standardDeviation*scalingFactor

    if standardDeviation_temp < 1e-14:
        standardDeviation_temp = 1e-14
    elif standardDeviation_temp > 0.1:
        standardDeviation_temp = 0.1
    
    return standardDeviation_temp

def getTemperature(temperature_minimum,temperature_maximum,step,step_maximum):
    
    return temperature_minimum * pow(temperature_maximum/temperature_minimum,(float(step)-1)/(float(step_maximum)-1))
    




