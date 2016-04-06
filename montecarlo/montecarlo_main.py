import copy
import random
import numpy as np
from printing import *
from tabulate import tabulate
from cost_functions import *
import time

class montecarlo:

    def __init__(self, _ham_num_, _zeta_, _momenta_,temperature_minimum,temperature_maximum,step_maximum,standard_deviation_maximum):

        self._ham_num_ = copy.deepcopy(_ham_num_)
        self._zeta_ = copy.deepcopy(_zeta_)
        self.momentum_values = _momenta_.getListOfMomenta()
        self.temperature_minimum = temperature_minimum
        self.temperature_maximum = temperature_maximum
        self.step_maximum = step_maximum
        self.nAccepted = 1.0
        self.nAttempts = 1.0
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
        self.berry = 0.0
        self.gap = 0.0

    def updateTemperature(self):

        self.temperature = self.temperature_minimum * pow(self.temperature_maximum/self.temperature_minimum,(float(self.step)-1)/(float(self.step_maximum)-1)) 

    def updateStandardDeviation(self):

        freq_accepted_ideal = 0.5
        scaling_factor = abs(self.freqAccepted - freq_accepted_ideal)
        if self.freqAccepted > freq_accepted_ideal: 
            self.standard_deviation = self.standard_deviation/scaling_factor
        else: 
            self.standard_deviation = self.standard_deviation*scaling_factor
        if self.standard_deviation < 1e-14: 
            self.standard_deviation = 1e-14
        elif self.standard_deviation > self.standard_deviation_maximum: 
            self.standard_deviation = self.standard_deviation_maximum

    def updateInternals(self,accepted_flag):

        self.nAttempts += 1
        if accepted_flag: 
            self.nAccepted += 1; 
            self.cost_function_current = self.cost_function_proposed; 
            self._zeta_.acceptMove()
        if self.cost_function_current < self.cost_function_minimum: 
            self.cost_function_minimum = self.cost_function_current
            self._zeta_minimum = copy.deepcopy(self._zeta_)
        self.freqAccepted = self.nAccepted/self.nAttempts
        self.step += 1

    def doSearch(self):
    
        stdscr = initialiseDisplay() 
        self.getCostFunction(True) 
        test = self.cost_function_current 
        f = open('berry.txt','w')
        while  self.step < self.step_maximum:    
            accepted_flag = False
            self.updateTemperature()
            self.updateStandardDeviation()
            self._zeta_.generateMove(self.standard_deviation)
            self.getCostFunction()
            self.updateAcceptanceProbability()
            if np.log(random.random()) <= self.acceptance_probability: accepted_flag = True
            self.updateInternals(accepted_flag)
            self.printState(stdscr,self._ham_num_.list_of_coefficients)
            f.write('%f %f\n' % (self.berry,self.gap))
        self.printState(stdscr,self._ham_num_.list_of_coefficients,True)
        f.close()
        killDisplay()

    def updateAcceptanceProbability(self):

        self.acceptance_probability = -self.temperature * (self.cost_function_proposed - self.cost_function_current)

    def getCostFunction(self,initial_flag=False):

        if initial_flag: 
            coefficient_values = self._zeta_.returnCurrentValues()
        else: 
            coefficient_values = self._zeta_.returnProposedValues()
        
        self.gap = pow(10,5)
        eSys = []

        for index,momentum_value in enumerate(self.momentum_values):
            eSys.append(self._ham_num_.calculateEigenvalues(coefficient_values + momentum_value))
            if abs(eSys[index][0][len(eSys[index][0][:])/2]) < self.gap: self.gap = (eSys[index][0][len(eSys[index][0][:])/2])
        
        
        self.berry = calculateBerryPhase(eSys)

        if initial_flag:    self.cost_function_current = np.abs(1-self.berry)
        else:               self.cost_function_proposed = np.abs(1-self.berry)

            
    def printState(self,stdscr,list_of_coefficients,end_flag=False):
        
        tables = [[],[],[]] 
        tables[0].append(['step','%d/%d'%(int(self.step),int(self.step_maximum))])
        tables[0].append(['temperature', '{0:.2f}'.format(self.temperature,self.temperature_maximum)])
        tables[0].append(['cost function (current)','%f'%(self.cost_function_current)])
        tables[0].append(['const funcion (minimum)','%f'%(self.cost_function_minimum)]) 
        tables[0].append(['berry','%f'%(self.berry)])
         
        for coeff,coeffval in zip(list_of_coefficients,self._zeta_.coefficients):
            tables[1].append([str(coeff),coeffval.real,coeffval.imag])
                  
        for coeff,coeffval in zip(list_of_coefficients,self._zeta_minimum.coefficients):
            tables[2].append([str(coeff),coeffval.real,coeffval.imag])
       
        headers = [[' ',' '],\
                ['Coefficient','  Current  ','  Current  '],\
                ['Coefficient','  Minimum  ','  Minimum  ']]

        if not end_flag:
            
            positions = [(0,0),(12,0),(40,0)]
            stdscr.addstr(positions[0][0],positions[0][1],tabulate(tables[0],headers=headers[0],tablefmt='grid'))
            stdscr.refresh()

            #for table,headers,position in zip(tables,headers,positions):
            #    stdscr.addstr(position[0],position[1],tabulate(table,headers=headers,tablefmt='grid'))
            #    stdscr.refresh()
        else:
            f = open(time.strftime('%d-%m-%y')+' '+time.strftime('%H:%M')+'.txt','a')
            for table,headers in zip(tables,headers): 
                f.write(tabulate(table,headers=headers,tablefmt='grid'))
                f.write('\n\n')

            f.close()
            
