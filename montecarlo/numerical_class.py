from sympy import *
import numpy
from eigensystem_functions import *

class numericalhamiltonian:

    def __init__(self,hamiltonian,list_of_coefficients,list_of_momenta):

        self.symbolic_hamiltonian = hamiltonian
        self.list_of_coefficients = list_of_coefficients
        self.lost_of_momenta = list_of_momenta
        self.hamiltonian = lambdify(list_of_coefficients + list_of_momenta,hamiltonian,"numpy")

    def calculateEigensystem(self,list_of_values):
        
        evals, evec = numpy.linalg.eig(self.hamiltonian(*list_of_values))
        evals_sort, evec_sort = sortEigensystem(evals,evec)

        return evals_sort, evec_sort
        
