import sys
#sys.path.insert(1,'../eigensystem/')
from eigensystem_functions import *
from numerical_class import *
import random
from sympy import *
import numpy as np
import math
import copy

def linspace(pMin,pMax,pInt):
    
    space = []
    space.append(pMin)
    inc = (pMax-pMin)/pInt
    for tik in range(0,pInt):
        space.append(space[tik]+inc)
    
    return space

class coeffcontainer:

    def __init__(self,real,imag,flag_real,flag_imag,real_free_flag,imag_free_flag):

        self.imag = imag
        self.real = real
        self.imag_flag = flag_imag
        self.real_flag = flag_real
        self.imag_free_flag = imag_free_flag
        self.real_free_flag = real_free_flag
        self.value = self.real + 1j * self.imag

    def assignValue(self,real,imag):

        self.real = real
        self.imag = imag
        self.value = self.real + 1j * self.imag

    def updateValue(self):
        
        self.value = self.real + 1j * self.imag

def getUniformSample(flag,width):

    if flag:
        return random.uniform(-width,width)    
    else:
        return 0.0
   
def getGaussianSample(flag,free_flag,standard_deviation,current_value):
    
    if free_flag:
        if flag: 
            return random.gauss(current_value,standard_deviation)
        else:
            return 0.0
    else: 
        return 1.0

class zeta:

    def __init__(self,zeta_input):

        self.coefficients = []
        self.coefficients_move = []
        
        for tik  in range(len(zeta_input)):
            flag_real = zeta_input[tik][0]
            flag_imag = zeta_input[tik][1]
            real_free_flag = zeta_input[tik][2]
            imag_free_flag = zeta_input[tik][3]
            self.coefficients.append(coeffcontainer(getUniformSample(flag_real,10.0),getUniformSample(flag_imag,10.0),flag_real,flag_imag,real_free_flag,imag_free_flag))
            self.coefficients_move.append(coeffcontainer(getUniformSample(flag_real,10.0),getUniformSample(flag_imag,10.0),flag_real,flag_imag,real_free_flag,imag_free_flag))
    
    def returnCurrentValues(self):

        return [ coeff.value for coeff in self.coefficients ]

    def returnProposedValues(self):

        return [ coeff.value for coeff in self.coefficients_move ]

    def findRedundant(self,_ham_num_,list_of_coefficients,list_of_momenta):
        
        _momentum_example_ = [ 1 for tik in range(len(list_of_momenta)) ]       
        ham_input_base = self.returnCurrentValues() + _momentum_example_
        base_evals = _ham_num_.calculateEigenvalues(ham_input_base)    
        
        shift = 1.123

        for tik in range(len(list_of_coefficients)):
            ham_input_shift = copy.copy(ham_input_base)
            ham_input_shift[tik] = ham_input_shift[tik] + shift
            shift_evals = _ham_num_.calculateEigenvalues(ham_input_shift) 
            if sum([abs(a-b) for a,b in zip(base_evals,shift_evals)]) == 0: self.coefficients[tik].real_free_flag = 0
            ham_input_shift = copy.copy(ham_input_base)
            ham_input_shift[tik] = ham_input_shift[tik] + 1j * shift
            shift_evals = _ham_num_.calculateEigenvalues(ham_input_shift) 
            if sum([abs(a-b) for a,b in zip(base_evals,shift_evals)]) == 0: self.coefficients[tik].imag_free_flag = 0 

    def generateMove(self,standard_deviation):
        
        count = 0

        for coeff in self.coefficients_move: 
            
            real = getGaussianSample(coeff.real_flag,coeff.real_free_flag,standard_deviation,coeff.real)
            imag = getGaussianSample(coeff.imag_flag,coeff.imag_free_flag,standard_deviation,coeff.imag)
            if count == 6: print(real,coeff.real,coeff.real_free_flag)
            coeff.assignValue(real,imag)
            count += 1
            
    def acceptMove(self):

        self.coefficients = copy.deepcopy(self.coefficients_move)

class momentum:

    def __init__(self,discretisation,spatial_dimension):
        
        pi = 3.14159
        self.momentum = [0 for tik in range(spatial_dimension)]
        self.space = linspace(0,pi,discretisation)
        self.dimension = spatial_dimension
        self.discretisation = discretisation

    def getRandom(self):

        return self.momentum

    def getListOfMomenta(self):

        if self.dimension == 1: return [[value] for value in self.space]
        elif self.dimension == 2: return [[value1,value2] for value1 in self.space for value2 in self.space ]
        elif self.dimension == 3: return [[value1,value2,value3] for value1 in self.space for value2 in self.space for value3 in self.space ]

