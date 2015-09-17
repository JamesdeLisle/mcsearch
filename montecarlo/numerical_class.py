from sympy import *
import scipy.linalg as LA
from eigensystem_functions import *

class numericalhamiltonian:

    def __init__(self,hamiltonian,list_of_coefficients,list_of_momenta):

        self.symbolic_hamiltonian = hamiltonian
        self.list_of_coefficients = list_of_coefficients
        self.lost_of_momenta = list_of_momenta
        self.hamiltonian = lambdify(list_of_coefficients + list_of_momenta,hamiltonian,"numpy")

    def calculateEigenvalues(self,list_of_values):
        
        evals = LA.eigvals(self.hamiltonian(*list_of_values))
        evals_sort = sortEigensystem(evals)

        return evals_sort
        
