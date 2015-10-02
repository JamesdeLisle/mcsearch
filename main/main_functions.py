from hamiltonian import *
from hamiltonian_functions import *
from general_functions import *
from symmetry_tests import  *
from sympy.physics.quantum import TensorProduct
from zeta_class import *

def generateHamiltonian(_param_): 
       
    vector_set = getListOfVectors(_param_.spatial_dimension,_param_.max_interaction_length)
    vector_symbols = getVectorSymbols(vector_set)
    list_of_species = getSpeciesSymbols(_param_.number_of_species)
    list_of_momenta = getMomentumSymbols(_param_.spatial_dimension) 
    list_of_terms,list_of_coefficients = getListOfTerms(vector_set,vector_symbols,list_of_species,_param_.superconducting_flag)
    H = hamiltonian(list_of_terms,list_of_coefficients,list_of_species)
    
    output = H.returnMatrix()

    return output,list_of_coefficients,list_of_momenta

def constrainHamiltonian(hamiltonian,_param_):
    
    if _param_.superconducting_flag:
        PH_matrix = TensorProduct(getPauli(1),eye(_param_.number_of_species))
        TR_matrix = TensorProduct(eye(_param_.number_of_species),getPauli(2))
    else:
        #TODO implement general construction of symmetry operators
        pass
    
    H_PH_inv = constructPHInvH(PH_matrix,hamiltonian)
    H_PH_TR_inv = constructTRInvH(TR_matrix,H_PH_inv)
     
    return H_PH_TR_inv

def generateZeta(hamiltonian,list_of_coefficients,list_of_momenta,_param_):

    zeta_input = [] 
    for variable in list_of_coefficients:
        if str(variable)[0] == 'm': zeta_input.append([1,0,1,1])
        else: zeta_input.append([1,1,1,1])
    
    _momenta_ = momentum(_param_.momentum_discretisation,_param_.spatial_dimension)
    _zeta_ = zeta(zeta_input)
    _zeta_.findRedundant(hamiltonian,list_of_coefficients,list_of_momenta)
    
    return _zeta_, _momenta_

