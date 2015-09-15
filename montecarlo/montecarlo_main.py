
import copy
import random
import numpy as np

class montecarlo:

    def __init__(self, _ham_num_, _zeta_, _momenta_,temperature_minimum,temperature_maximum,step_maximum,standard_deviation_maximum):

        self._ham_num_ = copy.deepcopy(_ham_num_)
        self._zeta_ = copy.deepcopy(_zeta_)
        self._momenta_ = copy.deepcopy(_momenta_)
        self.temperature_minimum = temperature_minimum
        self.temperature_maximum = temperature_maximum
        self.step_maximum = step_maximum
        self.nAccepted = 0.0
        self.nAttempts = 0.0
        self.freqAccepted = 0.0
        self.step = 1
        self.temperature = 0.0
        self.standard_deviation = 1.0
        self.standard_deviation_maximum = 0.1
        self.cost_function_current = 0.0
        self.cost_function_proposed = 0.0
        self.cost_function_minimum = pow(10,4)
        self._zeta_minimum = copy.deepcopy(self._zeta_)
        self.acceptance_probability = 0.0

    def updateTemperature(self):

        self.temperature = self.temperature_minimum * pow(self.temperature_maximum/self.temperature_minimum,(float(self.step)-1)/(float(self.step_maximum)-1)) 

    def updateStandardDeviation(self):

        freq_accepted_ideal = 0.5
        scaling_factor = abs(self.freqAccepted - freq_accepted_ideal)
        if self.freqAccepted > freq_accepted_ideal: self.standard_deviation = self.standard_deviation/scaling_factor
        else: self.standard_deviation = self.standard_deviation*scaling_factor
        if self.standard_deviation < 1e-14: self.standard_deviation = 1e-14
        elif self.standard_deviation > self.standard_deviation_maximum: self.standard_deviation = self.standard_deviation_maximum

    def updateInternals(self,accepted_flag):

        self.nAttempts += 1
        if accepted_flag: self.nAccepted += 1; self.cost_function_current = self.cost_function_proposed; self._zeta_.acceptMove()
        if self.cost_function_current < self.cost_function_minimum: self.cost_function_minimum = self.cost_function_current
        self.freqAccepted = self.nAccepted/self.nAttempts
        self.step += 1

    def doSearch(self):

        self.getCostFunction(True) 
         
        while  self.step < self.step_maximum:
            
            accepted_flag = False
            self.updateTemperature()
            self.updateStandardDeviation()
            self._zeta_.generateMove(self.standard_deviation)
            self.getCostFunction()
            self.updateAcceptanceProbability()
            if np.log(random.random()) <= self.acceptance_probability: accepted_flag = True
            self.updateInternals(accepted_flag)
            print(self.step)


    def updateAcceptanceProbability(self):

        self.acceptance_probability = -self.temperature * (self.cost_function_proposed - self.cost_function_current)

    def getCostFunction(self,initial_flag=False):

        momentum_values = self._momenta_.getListOfMomenta()
        if initial_flag: coefficient_values = self._zeta_.returnCurrentValues()
        else: coefficient_values = self._zeta_.returnProposedValues()
        
        for momentum_value in momentum_values:
            evals, evecs = self._ham_num_.calculateEigensystem(coefficient_values + momentum_value)

        if initial_flag: self.cost_function_current = 1
        else: self.cost_function_proposed = 1

            
