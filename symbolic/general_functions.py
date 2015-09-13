from sympy import *
from quadratic_operators import *

def getMomentumSymbols(spatial_dimension):

    output = []
    for tik in range(spatial_dimension):
        output.append(Symbol('p%d'%(tik),real=True))
    
    return output

def onlyHalfOfVec(vectorSet):
    
    vectorSetTemp = []
    zeroFlag = False
    
    for tik in range(0,len(vectorSet)):
        if vectorSet[tik] == zeros(len(vectorSet[tik]),1):
            zeroFlag = True
        if zeroFlag:
            vectorSetTemp.append(vectorSet[tik])
    
    return vectorSetTemp

def getListOfVectors(dimension,max_interaction_length):

    vector_set = []
    
    if dimension == 1:       
        for tik1 in range(-1,max_interaction_length+1):
            vector_set.append(Matrix([tik1]))

    elif dimension == 2:
        for tik1 in range(-1,max_interaction_length+1):
            for tik2 in range(-1,max_interaction_length+1):
                vector_set.append(Matrix([tik1,tik2]))

    elif dimension == 3:
        for tik1 in range(-1,max_interaction_length+1):
            for tik2 in range(-1,max_interaction_length+1):
                for tik3 in range(-1,max_interaction_length+1):
                    vector_set.append(Matrix([tik1,tik2,tik3]))
                                                        
    vector_set = onlyHalfOfVec(vector_set)
   
    return vector_set

def getVectorSymbols(vector_set):
    
    return symbols('s0:%d' % (len(vector_set)))

def getSpeciesSymbols(number_of_species):

    a,b,c,d = symbols('a b c d',real=True)
    types = [a,b,c,d]
    
    out = []
    for tik in range(number_of_species):
        out.append(types[tik])
            
    return out

def getListOfTerms(vector_set,vector_symbols,list_of_species,superconducting_flag):
    
    list_of_terms = []
    list_of_coefficients = []

    tSub = 0
    muSub = 0
    deltaSub = 0

    for tik in range(len(list_of_species)):
        list_of_terms.append(quadOp(Symbol('mu%d'%(muSub),real=True),list_of_species[tik],vector_set[0],vector_symbols[0],True,list_of_species[tik],vector_set[0],vector_symbols[0],False))
        list_of_coefficients.append(Symbol('mu%d'%(muSub),real=True))
        muSub = muSub + 1

    for tik1 in range(len(list_of_species)-1):
        list_of_terms.append(quadOp(Symbol('t%d'%(tSub)),list_of_species[tik1],vector_set[0],vector_symbols[0],True,list_of_species[tik1+1],vector_set[0],vector_symbols[0],False))
        list_of_coefficients.append(	Symbol('t%d'%(tSub)))
        tSub = tSub + 1
            

    if superconducting_flag:
        for tik1 in range(len(list_of_species)-1):
            list_of_terms.append(quadOp(Symbol('Delta%d'%(deltaSub)),list_of_species[tik1],vector_set[0],vector_symbols[0],False,list_of_species[tik1+1],vector_set[0],vector_symbols[0],False))
            list_of_coefficients.append(Symbol('Delta%d'%(deltaSub)))
            deltaSub = deltaSub + 1

    for tik in range(1,len(vector_set)):
        for tik1 in range(len(list_of_species)):
            for tik2 in range(len(list_of_species)):
                if True:
                    list_of_terms.append(quadOp(Symbol('t%d'%(tSub)),list_of_species[tik1],vector_set[0],vector_symbols[0],True,list_of_species[tik2],vector_set[tik],vector_symbols[tik],False))        
                    list_of_coefficients.append(Symbol('t%d'%(tSub)))
                    tSub = tSub + 1 	

                    if superconducting_flag:
                        list_of_terms.append(quadOp(Symbol('Delta%d'%(deltaSub)),list_of_species[tik1],vector_set[0],vector_symbols[0],False,list_of_species[tik2],vector_set[tik],vector_symbols[tik],False))
                        list_of_coefficients.append(Symbol('Delta%d'%(deltaSub)))
                        deltaSub = deltaSub + 1
 
    return list_of_terms,list_of_coefficients 

