import sys
sys.path.insert(1,'symbolic/')
sys.path.insert(1,'montecarlo/')
from hamiltonian import *
from hamiltonian_functions import *
from general_functions import *
from symmetry_tests import  *
from sympy.physics.quantum import TensorProduct
from zeta_class import *

def generateHamiltonian(dimension,max_interaction_length,number_of_species,superconducting_flag): 
       
    vector_set = getListOfVectors(dimension,max_interaction_length)
    vector_symbols = getVectorSymbols(vector_set)
    list_of_species = getSpeciesSymbols(number_of_species)
    list_of_momenta = getMomentumSymbols(dimension) 
    display(list_of_species)
    
    list_of_terms,list_of_coefficients = getListOfTerms(vector_set,vector_symbols,list_of_species,superconducting_flag)
    H = hamiltonian(list_of_terms,list_of_coefficients,list_of_species)
    #H.printReal()
    #H.printMatrix()
    #H.printMom()
    
    output = ImmutableMatrix(H.returnMatrix())
    output = output.rewrite(cos)
    output = output.rewrite(sin)
    #display(simplify(output))
    display(testPH(TensorProduct(getPauli(1),eye(number_of_species)),output))

    return output,list_of_coefficients,list_of_momenta

def constrainHamiltonian(hamiltonian,number_of_species,superconducting_flag):
    
    if superconducting_flag:
        #PH_matrix = TensorProduct(eye(number_of_species),getPauli(1))
        PH_matrix = TensorProduct(getPauli(1),eye(number_of_species))
        #TR_matrix = TensorProduct(getPauli(3),eye(number_of_species))
        TR_matrix = TensorProduct(eye(number_of_species),getPauli(2))
    else:
        #TODO implement general construction of symmetry operators
        pass
    
    H_PH_inv = constructPHInvH(PH_matrix,hamiltonian)
    H_PH_TR_inv = constructTRInvH(TR_matrix,H_PH_inv)
     
    return H_PH_TR_inv

def generateZeta(spatial_dimension,list_of_coefficients,list_of_momenta,hamiltonian,momentum_discretisation):

    zeta_input = [] 
    for variable in list_of_coefficients:
        if str(variable)[0] == 'm':
            zeta_input.append([1,0,1,1])
        else:
            zeta_input.append([1,1,1,1])
    
    _momenta_ = momentum(momentum_discretisation,spatial_dimension)
    _zeta_ = zeta(zeta_input)
    _zeta_.findRedundant(hamiltonian,list_of_coefficients,list_of_momenta,_momenta_)


spatial_dimension = 1
max_interaction_length = 3
number_of_species = 2 
superconducting_flag = True
momentum_discretisation = 30

H, list_of_coefficients,list_of_momenta = generateHamiltonian(spatial_dimension,max_interaction_length,number_of_species,superconducting_flag)
constrainHamiltonian(H,number_of_species,superconducting_flag)
generateZeta(spatial_dimension,list_of_coefficients,list_of_momenta,H,momentum_discretisation)
