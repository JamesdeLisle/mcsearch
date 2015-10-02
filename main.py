import sys
sys.path.insert(1,'symbolic/')
sys.path.insert(1,'montecarlo/')
sys.path.insert(1,'eigensystem/')
sys.path.insert(1,'main/')
sys.path.insert(1,'print/')
from montecarlo_main import *
from parameters_class import *
from main_functions import *

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
    _H_inv_ = constrainHamiltonian(_H_,_param_)
    #----------------------------------------------------------------#


    #----------------------------------------------------------------#
    # generate numerical hamiltonian
    _ham_num_ = numericalhamiltonian(_H_inv_,list_of_coefficients,list_of_momenta)
    _zeta_, _momenta_ = generateZeta(_ham_num_,list_of_coefficients,list_of_momenta,_param_)
    #----------------------------------------------------------------#


    #----------------------------------------------------------------#
    # confugure and run montecarlo search
    temperature_minimum = 1
    temperature_maximum = pow(10,6)
    step_maximum = pow(10,4)
    standard_deviation_maximum = 0.1

    _mc_ = montecarlo(_ham_num_, _zeta_, _momenta_,temperature_minimum,temperature_maximum,step_maximum,standard_deviation_maximum)
    _mc_.doSearch()
    #----------------------------------------------------------------#

main()




