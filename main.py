import sys
sys.path.insert(1,'symbolic/')
sys.path.insert(1,'montecarlo/')
sys.path.insert(1,'eigensystem/')
sys.path.insert(1,'main/')
sys.path.insert(1,'print/')
from hamiltonian import *
from hamiltonian_functions import *
from general_functions import *
from symmetry_tests import  *
from sympy.physics.quantum import TensorProduct
from zeta_class import *
from parameters_class import *
from montecarlo_main import *

def generateHamiltonian(_param_): 
       
    vector_set = getListOfVectors(_param_.spatial_dimension,_param_.max_interaction_length)
    vector_symbols = getVectorSymbols(vector_set)
    list_of_species = getSpeciesSymbols(_param_.number_of_species)
    list_of_momenta = getMomentumSymbols(_param_.spatial_dimension) 
    display(list_of_species)
    
    list_of_terms,list_of_coefficients = getListOfTerms(vector_set,vector_symbols,list_of_species,_param_.superconducting_flag)
    H = hamiltonian(list_of_terms,list_of_coefficients,list_of_species)
    
    output = H.returnMatrix()

    return output,list_of_coefficients,list_of_momenta

def constrainHamiltonian(hamiltonian,_param_,superconducting_flag):
    
    if superconducting_flag:
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

def main():

    #----------------------------------------------------------------#
    # configure initialisation variables
    spatial_dimension = 1
    max_interaction_length = 1
    number_of_species = 2
    superconducting_flag = True
    momentum_discretisation = 30
    _param_ = parameters(spatial_dimension,max_interaction_length,number_of_species,superconducting_flag,momentum_discretisation)
    #----------------------------------------------------------------#


    #----------------------------------------------------------------#
    # generate and restrict symbolic hamiltonian
    _H_, list_of_coefficients,list_of_momenta = generateHamiltonian(_param_)
    _H_inv_ = constrainHamiltonian(_H_,_param_,superconducting_flag)
    #----------------------------------------------------------------#


    #----------------------------------------------------------------#
    # generate numerical hamiltonian
    _ham_num_ = numericalhamiltonian(_H_inv_,list_of_coefficients,list_of_momenta)
    _zeta_, _momenta_ = generateZeta(_ham_num_,list_of_coefficients,list_of_momenta,_param_)
    #----------------------------------------------------------------#


    #----------------------------------------------------------------#
    # confugure and run montecarlo search
    temperature_minimum = 1
    temperature_maximum = pow(10,5)
    step_maximum = pow(10,4)
    standard_deviation_maximum = 0.1

    _mc_ = montecarlo(_ham_num_, _zeta_, _momenta_,temperature_minimum,temperature_maximum,step_maximum,standard_deviation_maximum)
    _mc_.doSearch()
    #----------------------------------------------------------------#

main()




