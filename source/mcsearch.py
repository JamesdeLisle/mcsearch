from numpy import array, log
from printFunctions import *
from mcsearch_functions import *
from generate_hamiltonian import *
import random
from lambda_class import *
import time

print_flag = 1

def monteCarloSearch():
    
    if print_flag:
        stdscr = initialiseDisplay()
    
    #TRmatrix = array([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]])
    #TRmatrix = array([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
    TRmatrix = array([[1,0,0,0],[0,-1,0,0],[0,0,1,0],[0,0,0,-1]])
    
    coefficient_names = {'mu0':[1,0,1,1],\
            'mu1':[1,0,1,1],'t0':[1,1,1,1],\
            'Delta0':[1,1,1,1],'t1':[1,1,1,1],\
            'Delta1':[1,1,1,1],'t2':[1,1,1,1],\
            'Delta2':[1,1,1,1],'t3':[1,1,1,1],\
            'Delta3':[1,1,1,1],'t4':[1,1,1,1],\
            'Delta4':[1,1,1,1]}

    momentum_discretisation = 50
    _Lambda_ = Lambda(coefficient_names) 
    _Lambda_.findFree(TRmatrix,momentum_discretisation)
    costFunction_current = getCostFunction(_Lambda_.coefficients,TRmatrix,momentum_discretisation)  
    costFunction_minimum = costFunction_current  
    Lambda_minimum = _Lambda_	
    nAccepted = 0.0                          
    nAttempts = 0.0                      
    freqAccepted = 0.0                          
    temperature = 0.0 
    temperature_minimum = 1
    temperature_maximum = 2 * pow(10,4)                     
    standard_deviation = 1.0              
    step = 1
    step_maximum = 5 * pow(10,4) 

    while step < step_maximum:
                                                                                                    
        #Lambda_current = renormaliseLambda(Lambda_current)
        temperature = getTemperature(temperature_minimum,temperature_maximum,step,step_maximum)
        nAttempts = nAttempts + 1.0 
        standard_deviation = updateStandardDeviation(freqAccepted,standard_deviation) 
        _Lambda_.generateMove(standard_deviation)
        #for coeff in _Lambda_.coefficients:
        #    print(_Lambda_.coefficients[coeff].imag_flag)
        #    print(_Lambda_.coefficients_move[coeff].imag_flag)
        #    print(' ')
        #print('----------')
        costFunction_move = getCostFunction(_Lambda_.coefficients_move,TRmatrix,momentum_discretisation)
        
        acceptanceProbability = -temperature*(costFunction_move-costFunction_current)
        draw = random.random()
        
        if print_flag:
            printLambda(stdscr,_Lambda_)
            printTable(stdscr, step,step_maximum,temperature,costFunction_current,costFunction_minimum)
        
        key = stdscr.getch()
        if key == ord('q'):
            killDisplay()

        if log(draw) <= acceptanceProbability: 
            #print('hi')
            _Lambda_.acceptMove()
            costFunction_current = costFunction_move
            nAccepted = nAccepted + 1.0
          
        freqAccepted = nAccepted/nAttempts               
        #print(costFunction_move) 
        #print(costFunction_current) 
        #print('-------------')

        if costFunction_current < costFunction_minimum:
            costFunction_minimum = costFunction_current
            Lambda_minimum = _Lambda_
 
        step = step + 1 
    
    killDisplay()
    #printLambdaFinal(Lambda_current)
    #printRunToFile(costFunction_minimum,Lambda_minimum,Lambda_current,costFunction_current, \
    #        TRmatrix,temperature_maximum,momentum_discretisation)
#####################################################
monteCarloSearch()

