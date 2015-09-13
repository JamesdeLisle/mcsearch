import random
from mcsearch_functions import *

class coeffContainer:

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

class Lambda:

    def __init__(self,coefficient_names):

        self.coefficients = {}
        self.coefficients_move = {}
        for coeff in coefficient_names:
            flag_real = coefficient_names[coeff][0]
            flag_imag = coefficient_names[coeff][1]
            real_free_flag = coefficient_names[coeff][2]
            imag_free_flag = coefficient_names[coeff][3]
            self.coefficients[coeff] = coeffContainer(getUniformSample(flag_real,10.0),getUniformSample(flag_imag,10.0),flag_real,flag_imag,real_free_flag,imag_free_flag)
            self.coefficients_move[coeff] = coeffContainer(getUniformSample(flag_real,10.0),getUniformSample(flag_imag,10.0),flag_real,flag_imag,real_free_flag,imag_free_flag)
    
    def findFree(self,TRmatrix,momentum_discretisation):
 
        costFunction_compare = getCostFunction(self.coefficients,TRmatrix,momentum_discretisation)  
        costFunction_test = 0.0              
        shift = 5
         
        for coeff in self.coefficients:
            self.coefficients[coeff].real = self.coefficients[coeff].real + shift
            self.coefficients[coeff].updateValue() 
            costFunction_test = getCostFunction(self.coefficients,TRmatrix,momentum_discretisation)
            
            if costFunction_test - costFunction_compare == 0:
                
                self.coefficients[coeff].real_free_flag = 0
                self.coefficients_move[coeff].real_free_flag = 0
                self.coefficients[coeff].real = self.coefficients[coeff].real - shift
            else:
                self.coefficients[coeff].real = self.coefficients[coeff].real - shift
            
            self.coefficients[coeff].updateValue() 
            self.coefficients[coeff].imag = self.coefficients[coeff].imag + shift
            self.coefficients[coeff].updateValue()
            costFunction_test = getCostFunction(self.coefficients,TRmatrix,momentum_discretisation)
            
            if costFunction_test - costFunction_compare == 0:
                self.coefficients[coeff].imag_flag = 0
                self.coefficients_move[coeff].imag_flag = 0
                self.coefficients[coeff].imag = self.coefficients[coeff].imag - shift
            else:
                self.coefficients[coeff].imag = self.coefficients[coeff].imag - shift
            
            self.coefficients[coeff].updateValue()

    def generateMove(self,standard_deviation):

        for coeff in self.coefficients_move: 
            real = getGaussianSample(self.coefficients[coeff].real_flag,self.coefficients[coeff].real_free_flag,standard_deviation,self.coefficients[coeff].real)
            imag = getGaussianSample(self.coefficients[coeff].imag_flag,self.coefficients[coeff].imag_free_flag,standard_deviation,self.coefficients[coeff].imag)
            self.coefficients_move[coeff].assignValue(real,imag)
            


    def acceptMove(self):

        self.coefficients = self.coefficients_move




        


