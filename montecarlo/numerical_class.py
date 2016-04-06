from sympy import *
import scipy.linalg as LA
import numpy as np
from eigensystem_functions import *

class numericalhamiltonian:

    def __init__(self,hamiltonian,list_of_coefficients,list_of_momenta):

        self.symbolic_hamiltonian = hamiltonian
        self.list_of_coefficients = list_of_coefficients
        self.lost_of_momenta = list_of_momenta
        self.hamiltonian = lambdify(list_of_coefficients + list_of_momenta,hamiltonian,"numpy")

    def calculateEigenvalues(self,list_of_values):
        
        evals, evecs = LA.eig(self.hamiltonian(*list_of_values))
        evals = [ np.real(entry).item(0) for entry in evals]
        track = [ i for i in range(0,len(evals))]
        evals_sorted = sorted(evals)
        evecs_sorted = []

        for index,entry in enumerate(evals):
            for sorted_index,sorted_entry in enumerate(evals_sorted):
                if entry == sorted_entry:
                    track[index] = sorted_index
        
        for entry in track:
            evecs_sorted.append(evecs[entry])

        return evals_sorted, evecs_sorted
        
