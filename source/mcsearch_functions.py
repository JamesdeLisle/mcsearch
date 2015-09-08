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

def getCostFunction(Lambda,TRMatrix,momentum_discretisation):
    
    pi = 3.14159
    momentumSpace = linspace(0,pi,momentum_discretisation)
    costFunction = 0.0      

    for tik1 in range(0,len(momentumSpace)):
        Ham_PlusP = getHam(Lambda,momentumSpace[tik1])
        Ham_MinusP = getHam(Lambda,-momentumSpace[tik1])
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

def getActiveFlags(Lambda_current,TRmatrix,momentum_discretisation):
    
    activeFlags = [],[] 
    Lambda_test = Lambda_current 
    costFunction_compare = getCostFunction(Lambda_current,TRmatrix,momentum_discretisation)  
    costFunction_test = 0.0              
    shift = 2.51            
    
    for tik1 in range(0,12):
        for tik2 in range(0,2):
            Lambda_test[tik2][tik1] = Lambda_test[tik2][tik1] + shift
            costFunction_test = getCostFunction(Lambda_test,TRmatrix,momentum_discretisation)

            if absolute(costFunction_test - costFunction_compare) > 0.0:
                activeFlags[tik2].append(True)
            else:
                activeFlags[tik2].append(False)
                                                                                                                                                    
            Lambda_test[tik2][tik1] = Lambda_test[tik2][tik1] - shift

    return activeFlags

def getLambda(Lambda_current,standard_deviation,activeFlags):

    Lambda_move = [],[]

    for tik1 in range(0,12):
        for tik2 in range(0,2):
            if activeFlags[tik2][tik1]:
                if tik1 < 2 and tik2 == 1:
                    Lambda_move[tik2].append(0.0)	
                else:
                    Lambda_move[tik2].append(random.gauss(Lambda_current[tik2][tik1],standard_deviation))
            else:
                if tik1 < 2 and tik2 == 1:
                    Lambda_move[tik2].append(0.0)	
                else:
                    Lambda_move[tik2].append(1.0)

    return Lambda_move

def renormaliseLambda(Lambda):
   
    sum_real = sum(Lambda[0][:])
    sum_imag = sum(Lambda[1][:])
 
    for tik1 in range(0,2):
        for tik2 in range(0,12):
            if tik1 == 0:
                Lambda[tik1][tik2] = Lambda[tik1][tik2]/sum_real
            else:
                Lambda[tik1][tik2] = Lambda[tik1][tik2]/sum_imag
 
    return Lambda

def getTemperature(temperature_minimum,temperature_maximum,step,step_maximum):
    
    return temperature_minimum * pow(temperature_maximum/temperature_minimum,(float(step)-1)/(float(step_maximum)-1))
    




