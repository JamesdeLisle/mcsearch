from __future__ import print_function
import random
from numpy import real, dot, absolute, array, conjugate, linalg, set_printoptions, sqrt, trace, divide
from math import cos,sin,exp, pow
import cmath as cm
import matplotlib.pyplot as plt
from printFunctions import *
import sys
import time

def linspace(pMin,pMax,pInt):
    
    space = []
    space.append(pMin)
    inc = (pMax-pMin)/pInt
    for tik in range(0,pInt):
        space.append(space[tik]+inc)
    
    return space

def getHam(Lambda,p):

    p0 = p
    mu0	    = Lambda[0][0] + 1j*Lambda[1][0]
    mu1	    = Lambda[0][1] + 1j*Lambda[1][1] 
    t0	    = Lambda[0][2] + 1j*Lambda[1][2]
    Delta0  = Lambda[0][3] + 1j*Lambda[1][3]
    t1	    = Lambda[0][4] + 1j*Lambda[1][4]
    Delta1  = Lambda[0][5] + 1j*Lambda[1][5]
    t2	    = Lambda[0][6] + 1j*Lambda[1][6]
    Delta2  = Lambda[0][7] + 1j*Lambda[1][7]
    t3	    = Lambda[0][8] + 1j*Lambda[1][8]
    Delta3  = Lambda[0][9] + 1j*Lambda[1][9]
    t4	    = Lambda[0][10] + 1j*Lambda[1][10]
    Delta4  = Lambda[0][11] + 1j*Lambda[1][11] 
    I = 1j
    
    Ham = array([[mu0 + t1*cm.exp(I*p0)/2 + cm.exp(-I*p0)*conjugate(t1)/2, t0/2 + t2*cm.exp(I*p0)/2 + cm.exp(-I*p0)*conjugate(t3)/2, -Delta1*cm.exp(I*p0)/2 + cm.exp(-I*p0)*conjugate(Delta1)/2, -Delta0/2 - Delta2*cm.exp(I*p0)/2 + cm.exp(-I*p0)*conjugate(Delta3)/2], [t3*cm.exp(I*p0)/2 + conjugate(t0)/2 + cm.exp(-I*p0)*conjugate(t2)/2, mu1 + t4*cm.exp(I*p0)/2 + cm.exp(-I*p0)*conjugate(t4)/2, -Delta3*cm.exp(I*p0)/2 + conjugate(Delta0)/2 + cm.exp(-I*p0)*conjugate(Delta2)/2, -Delta4*cm.exp(I*p0)/2 + cm.exp(-I*p0)*conjugate(Delta4)/2], [-Delta1*cm.exp(-I*p0)/2 + cm.exp(I*p0)*conjugate(Delta1)/2, -Delta3*cm.exp(-I*p0)/2 + cm.exp(I*p0)*conjugate(Delta2)/2 + conjugate(Delta0)/2, -mu0 - t1*cm.exp(I*p0)/2 - cm.exp(-I*p0)*conjugate(t1)/2, -t0/2 - t2*cm.exp(I*p0)/2 - cm.exp(-I*p0)*conjugate(t3)/2], [-Delta0/2 - Delta2*cm.exp(-I*p0)/2 + cm.exp(I*p0)*conjugate(Delta3)/2, -Delta4*cm.exp(-I*p0)/2 + cm.exp(I*p0)*conjugate(Delta4)/2, -t3*cm.exp(I*p0)/2 - conjugate(t0)/2 - cm.exp(-I*p0)*conjugate(t2)/2, -mu1 - t4*cm.exp(I*p0)/2 - cm.exp(-I*p0)*conjugate(t4)/2]])

    #Ham = array([[mu0 + t1*cos(p0), t0/2.0 + t2*cm.exp(I*p0)/2.0 + t3*cm.exp(-I*p0)/2.0,-I*Delta1*sin(p0), -Delta0/2.0 - Delta2*cm.exp(I*p0)/2 + Delta3*cm.exp(-I*p0)/2.0],[t0/2.0 + t2*cm.exp(-I*p0)/2.0 + t3*cm.exp(I*p0)/2.0, mu1 + t4*cos(p0), Delta0/2.0 + Delta2*cm.exp(-I*p0)/2.0 - Delta3*cm.exp(I*p0)/2.0, -I*Delta4*sin(p0)], [I*Delta1*sin(p0),Delta0/2.0 + Delta2*cm.exp(I*p0)/2.0 - Delta3*cm.exp(-I*p0)/2.0, -mu0 - t1*cos(p0),-t0/2.0 - t2*cm.exp(I*p0)/2 - t3*cm.exp(-I*p0)/2.0], [-Delta0/2.0 - Delta2*cm.exp(-I*p0)/2.0 + Delta3*cm.exp(I*p0)/2.0, I*Delta4*sin(p0), -t0/2.0 - t2*cm.exp(-I*p0)/2.0 - t3*cm.exp(I*p0)/2.0,-mu1 - t4*cos(p0)]])
    
    return Ham

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
    #print(Lambda) 
    for tik1 in range(0,2):
        for tik2 in range(0,12):
            if tik1 == 0:
                Lambda[tik1][tik2] = Lambda[tik1][tik2]/sum_real
            else:
                Lambda[tik1][tik2] = Lambda[tik1][tik2]/sum_imag
    #print(Lambda)
    return Lambda

def monteCarloSearch():
    
    #TRmatrix = array([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]])
    TRmatrix = array([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])

    momentum_discretisation = 50
    Lambda_current = [],[]
    
    for tik in range(0,12):
        if tik < 2:
            Lambda_current[0].append(random.uniform(-10.0,10.0))
            Lambda_current[1].append(0.0)
        else:
            Lambda_current[0].append(random.uniform(-10.0,10.0))
            Lambda_current[1].append(random.uniform(-10.0,10.0))

    for tik in range(0,12):
        print('%f + %fi'%(Lambda_current[0][tik],Lambda_current[1][tik]))

    activeFlags = getActiveFlags(Lambda_current,TRmatrix,momentum_discretisation)

    for tik in range(0,12):
        print('%s, %s'%(activeFlags[0][tik],activeFlags[1][tik]))

    costFunction_current = getCostFunction(Lambda_current,TRmatrix,momentum_discretisation)  
    costFunction_minimum = costFunction_current	    
    Lambda_minimum = Lambda_current	
    nAccepted = 0.0                          
    nAttempts = 0.0                      
    freqAccepted = 0.0                          
    temperature = 0.0 
    temperature_minimum = 1
    temperature_maximum = 200                         
    standardDeviation = 1.0              
    step = 1
    step_maximum = pow(10,6) 
    
    while step < step_maximum:
                                                                                                    
        #Lambda_current = renormaliseLambda(Lambda_current)
        temperature = temperature_minimum * pow(temperature_maximum/temperature_minimum,(float(step)-1)/(float(step_maximum)-1))
        nAttempts = nAttempts + 1.0 
        standardDeviation = updateStandardDeviation(freqAccepted,standardDeviation) 
        Lambda_move = getLambda(Lambda_current,standardDeviation,activeFlags) 
        costFunction_move = getCostFunction(Lambda_move,TRmatrix,momentum_discretisation) 
        acceptanceProbability = min(1,exp(-temperature*(costFunction_move-costFunction_current))) 
        draw = random.random()
        #print(Lambda_move) 
        sys.stdout.write('\r' + '%d/%d || %f || %f || %f || %f    ' %(step,step_maximum,temperature,costFunction_move,costFunction_current,acceptanceProbability))
        sys.stdout.flush()

        #print(' ')
        #print(costFunction_move)
        #print(costFunction_current) 
        #print(acceptanceProbability)

        if draw <= acceptanceProbability:   
            Lambda_current = Lambda_move
            costFunction_current = costFunction_move
            nAccepted = nAccepted + 1.0
                                                                                                                                        
        freqAccepted = nAccepted/nAttempts               
 
        if costFunction_current < costFunction_minimum:
            costFunction_minimum = costFunction_current
            Lambda_minimum = Lambda_current
 
        step = step + 1
        
    print('\n')
    printLambdaToConsole(Lambda_current)
    print(costFunction_current)
    printRunToFile(costFunction_minimum,Lambda_minimum,Lambda_current,costFunction_current, \
                    TRmatrix,temperature_maximum,temperature_increment,momentum_discretisation)
#####################################################
monteCarloSearch()

