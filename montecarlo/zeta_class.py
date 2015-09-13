import random
from sympy import *
import numpy as np

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

    def findRedundant(self,hamiltonian,list_of_coefficients,list_of_momenta,_momentum_):
        
        coeff = list_of_coefficients + list_of_momenta
        num_ham = lambdify(coeff,hamiltonian,"numpy")
        ham_input = self.returnCurrentValues() + _momentum_.getRandom()
        base_evals, base_evec = np.linalg.eig(num_ham(*ham_input))
        print(base_evals,base_evec)

    def generateMove(self,standard_deviation):

        for coeff in self.coefficients_move: 
            real = getGaussianSample(self.coefficients[coeff].real_flag,self.coefficients[coeff].real_free_flag,standard_deviation,self.coefficients[coeff].real)
            imag = getGaussianSample(self.coefficients[coeff].imag_flag,self.coefficients[coeff].imag_free_flag,standard_deviation,self.coefficients[coeff].imag)
            self.coefficients_move[coeff].assignValue(real,imag)
            
    def acceptMove(self):

        self.coefficients = self.coefficients_move

class momentum:

    def __init__(self,discretisation,spatial_dimension):
        
        pi = 3.14159
        self.momentum = [0 for tik in range(spatial_dimension)]
        self.space = linspace(0,pi,discretisation)
        self.dimension = spatial_dimension
        self.discretisation = discretisation

    def getRandom(self):
        
        for tik in range(self.dimension):
            self.momentum[tik] = self.space[random.randrange(0,self.discretisation-1)]

        return self.momentum
