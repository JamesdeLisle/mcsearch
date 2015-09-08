from numpy import array, log
from printFunctions import *
from mcsearch_functions import *
from generate_hamiltonian import *
import random

def monteCarloSearch():
    
    stdscr = initialiseDisplay()
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

    activeFlags = getActiveFlags(Lambda_current,TRmatrix,momentum_discretisation)
    printLambdaFlags(stdscr,activeFlags)
 
    costFunction_current = getCostFunction(Lambda_current,TRmatrix,momentum_discretisation)  
    costFunction_minimum = costFunction_current	    
    Lambda_minimum = Lambda_current	
    nAccepted = 0.0                          
    nAttempts = 0.0                      
    freqAccepted = 0.0                          
    temperature = 0.0 
    temperature_minimum = 1
    temperature_maximum = 2 * pow(10,5)                     
    standardDeviation = 1.0              
    step = 1
    step_maximum = 5 * pow(10,4) 

    while step < step_maximum:
                                                                                                    
        #Lambda_current = renormaliseLambda(Lambda_current)
        temperature = getTemperature(temperature_minimum,temperature_maximum,step,step_maximum)
        nAttempts = nAttempts + 1.0 
        standardDeviation = updateStandardDeviation(freqAccepted,standardDeviation) 
        Lambda_move = getLambda(Lambda_current,standardDeviation,activeFlags) 
        costFunction_move = getCostFunction(Lambda_move,TRmatrix,momentum_discretisation) 
        acceptanceProbability = -temperature*(costFunction_move-costFunction_current)
        draw = random.random()
        printLambda(stdscr,Lambda_current)
        printTable(stdscr, step,step_maximum,temperature,costFunction_current,costFunction_minimum)
        

        if log(draw) <= acceptanceProbability:   
            Lambda_current = Lambda_move
            costFunction_current = costFunction_move
            nAccepted = nAccepted + 1.0
                                                                                                                                        
        freqAccepted = nAccepted/nAttempts               
 
        if costFunction_current < costFunction_minimum:
            costFunction_minimum = costFunction_current
            Lambda_minimum = Lambda_current
 
        step = step + 1
    
    killDisplay()
    printLambdaFinal(Lambda_current)
    printRunToFile(costFunction_minimum,Lambda_minimum,Lambda_current,costFunction_current, \
            TRmatrix,temperature_maximum,momentum_discretisation)
#####################################################
monteCarloSearch()

