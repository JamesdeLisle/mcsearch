from sympy import  *
from IPython.display import display
import copy
from quadratic_operators import *
from hamiltonian_functions import *

init_printing()

class hamiltonian:

    # Initialise internal list of terms based on input list
    # and create fourier transformed version of list
    def __init__(self,listofterms,listOfCoeff,listOfSpecies):
            
        self.listOfSpecies = listOfSpecies
        self.listOfCoeff = listOfCoeff
        self.numOfSpecies = len(self.listOfSpecies)
        self.listOfTerms = listofterms
        self.listOfTerms_real = []
        
        superconducting_flag = checkForSuperconducting(self.listOfTerms)        
        self.matrixDim, self.symbolicMatrix = assignMatrix(superconducting_flag,self.numOfSpecies)
        self.listOfTerms_real = returnHermitianForm(self.listOfTerms) 
        self.listOfTerms_momentum = returnMomentumForm(self.listOfTerms_real,superconducting_flag)
        self.listOfTerms_momentum = orderCanonically(self.listOfTerms_momentum) 
        self.matrixForm = returnMatrixForm(self.matrixDim,self.numOfSpecies,self.listOfSpecies,superconducting_flag,self.listOfTerms_momentum)
        
        for tik1 in range(self.matrixDim):
            for tik2 in range(self.matrixDim):
                self.symbolicMatrix[tik1,tik2] = copy.deepcopy(self.matrixForm[tik1][tik2].coeff)

    def printReal(self):
            
        self.RealSymbolic = Symbol('1')
        
        for entry in range(0,len(self.listOfTerms_real)):
            self.RealSymbolic = self.RealSymbolic + self.listOfTerms_real[entry].symbol
        
        display(self.RealSymbolic-Symbol('1'))
    #################################################

    #################################################
    # Prints the momentum space Hamiltonian
    def printMom(self):
            
        self.MomSymbolic = Symbol('1')
    
        for entry in range(0,len(self.listOfTerms_momentum)):
            self.MomSymbolic = self.MomSymbolic + self.listOfTerms_momentum[entry].symbol
        
        display(self.MomSymbolic-Symbol('1'))
    #################################################

    #################################################
    # Prints the matrix form of the momentum Hamiltonian 
    def printMatrix(self): 
        
        display(self.symbolicMatrix)
    #################################################

    #################################################
    # Returns matrix of sympy symbols
    def returnMatrix(self):
        return self.symbolicMatrix
    #################################################

