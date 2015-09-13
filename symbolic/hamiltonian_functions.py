from sympy import *
from quadratic_operators import *
from IPython.display import display
import copy

def compare(a,b):
    return a == b

def checkForSuperconducting(listOfTerms):

    flag = False
    for term in listOfTerms:            
        if term.DAGboola and term.DAGboolb:
            flag = True
            break
        elif not term.DAGboola and not term.DAGboolb:
            flag = True
            break
        else:
            pass
    
    return flag

def assignMatrix(superconducting_flag,numOfSpecies):

    size = 0
    if superconducting_flag:
        size = 2*numOfSpecies 
    else:
        size = numOfSpecies

    return size,zeros(size,size)

def returnHermitianForm(listOfTerms):
    
    output = [] 
    for term in listOfTerms:
        output.append(copy.deepcopy(term))
        herm = copy.deepcopy(term)
        herm.HC()
        output.append(herm) 

    return output

def returnMomentumForm(listOfTerms_real,superconducting_flag):

    listOfTerms_momentum_plus = copy.deepcopy(listOfTerms_real)
    listOfTerms_momentum_minus = copy.deepcopy(listOfTerms_real) 
    map(lambda quadOp: quadOp.FT(), listOfTerms_momentum_plus)
    map(lambda quadOp: quadOp.FT(), listOfTerms_momentum_minus)
    map(lambda quadOp: quadOp.minusP(), listOfTerms_momentum_minus)

    if superconducting_flag:
        final = listOfTerms_momentum_plus + listOfTerms_momentum_minus
        for term in final:
            term.coeff = term.coeff/2
        return final
    else:
        return listOfTerms_momentum_plus
    
def orderCanonically(listOfTerms_momentum):
    
    #TODO write member method in quadOp for commutation

    for entry in listOfTerms_momentum: 
        if entry.indexaSym == entry.indexbSym: 
            if entry.indexaSym == -Symbol('p',real=True):
                if entry.DAGboola:
                    entry.DAGboola = not entry.DAGboola
                    entry.DAGboolb = not entry.DAGboolb
                    temptypea = entry.typea 
                    entry.typea = copy.deepcopy(entry.typeb)
                    entry.typeb = copy.deepcopy(temptypea)
                    entry.coeff = -entry.coeff
            else:
                if entry.DAGboolb:
                    entry.DAGboola = not entry.DAGboola
                    entry.DAGboolb = not entry.DAGboolb
                    temptypea = entry.typea 
                    entry.typea = copy.deepcopy(entry.typeb)
                    entry.typeb = copy.deepcopy(temptypea)
                    entry.coeff = -entry.coeff
        else:
            if entry.DAGboola:
                if entry.indexaSym == -Symbol('p',real=True):
                    entry.indexaSym = -entry.indexaSym
                    entry.indexbSym = -entry.indexbSym
                    temptypea = entry.typea 
                    entry.typea = copy.deepcopy(entry.typeb)
                    entry.typeb = copy.deepcopy(temptypea)
                    entry.coeff = -entry.coeff
            else:
                 if entry.indexaSym == Symbol('p',real=True):
                    entry.indexaSym = -entry.indexaSym
                    entry.indexbSym = -entry.indexbSym
                    temptypea = entry.typea 
                    entry.typea = copy.deepcopy(entry.typeb)
                    entry.typeb = copy.deepcopy(temptypea)
                    entry.coeff = -entry.coeff
        
        entry.setSymbol() 

    return listOfTerms_momentum

def returnMatrixForm(matrixDim,numOfSpecies,listOfSpecies,superconducting_flag,listOfTerms_momentum):
    
    #TODO rewite this function in a nicer way
    simp = []
    for tik in range(0,matrixDim):
        simp.append([])

    for tik1 in range(0,numOfSpecies):
        for tik2 in range(0,numOfSpecies):
            simp[tik1].append(quadOp(Symbol('1',real=True), \
                    listOfSpecies[tik1],[],Symbol('p',real=True), \
                    True,listOfSpecies[tik2],[],Symbol('p',real=True),False))
        if superconducting_flag:
            for tik2 in range(0,numOfSpecies):
                simp[tik1].append(quadOp(Symbol('1',real=True), \
                        listOfSpecies[tik1],[],Symbol('p',real=True), \
                        True,listOfSpecies[tik2],[],-Symbol('p',real=True),True))
    for tik1 in range(numOfSpecies,matrixDim):
            for tik2 in range(0,numOfSpecies):
                simp[tik1].append(quadOp(Symbol('1',real=True), \
                        listOfSpecies[tik1-numOfSpecies],[], \
                        -Symbol('p',real=True),False,listOfSpecies[tik2], \
                        [],Symbol('p',real=True),False))
            if superconducting_flag:
                for tik2 in range(0,numOfSpecies):
                    simp[tik1].append(quadOp(Symbol('1',real=True), \
                            listOfSpecies[tik1-numOfSpecies],[], \
                            -Symbol('p',real=True),False,listOfSpecies[tik2], \
                            [],-Symbol('p',real=True),True))
    
    for tik1 in range(0,matrixDim):
        for tik2 in range(0,matrixDim):             
            for tik3 in range(0,len(listOfTerms_momentum)):

                bool1 = compare(simp[tik1][tik2].DAGboola,listOfTerms_momentum[tik3].DAGboola)
                bool2 = compare(simp[tik1][tik2].DAGboolb,listOfTerms_momentum[tik3].DAGboolb)
                bool3 = compare(simp[tik1][tik2].indexaSym,listOfTerms_momentum[tik3].indexaSym)
                bool4 = compare(simp[tik1][tik2].indexbSym,listOfTerms_momentum[tik3].indexbSym)
                bool5 = compare(simp[tik1][tik2].typea,listOfTerms_momentum[tik3].typea)
                bool6 = compare(simp[tik1][tik2].typeb,listOfTerms_momentum[tik3].typeb)
                
                if bool1 * bool2 * bool3 * bool4 * bool5 * bool6:    
                    simp[tik1][tik2].coeff = simp[tik1][tik2].coeff + listOfTerms_momentum[tik3].coeff

    for tik1 in range(0,matrixDim):
        for tik2 in range(0,matrixDim):
            simp[tik1][tik2].coeff = simp[tik1][tik2].coeff - Symbol('1',real=True)

    return simp

def simplifyExp(F):

    fSy = [] 
    X = F
    for tik in range(0,len(X.args)):
        if len(X.args[tik].args) != 0:
            for tik1 in range(0,len(X.args[tik].args)):
                if type(X.args[tik].args[tik1]) == type(Symbol('a')):
                    fSy.append(X.args[tik].args[tik1])

    for tik in range(0,len(fSy)):
        if type(fSy[tik]) == type(Symbol('a')):
            X = collect(X,fSy[tik])

    f1Sy = []
    for tik in range(0,len(X.args)):
        if len(X.args[tik].args) != 0:
            f1Sy.append(X.args[tik].args[1])
        else:
            f1Sy.append(X.args[1])
               
    for tik in range(0,len(f1Sy)):
        if type(f1Sy[tik]) == type(Symbol('a')+Symbol('b')):
            X = X.subs(f1Sy[tik],f1Sy[tik].rewrite(cos)) 
    
    return X

def listGenerate(vectorSet, vectorLabels, listOfSpecies, isSuper):
	
    v = vectorSet
    s = vectorLabels
    los = listOfSpecies
    lot = []
    
    for tik in range(0,len(listOfSpecies)):
        lot.append(quadOp(Symbol('mu',real=True),los[tik],v[0],s[0],True,los[tik],v[0],s[0],False))
       
    for tik1 in range(0,len(listOfSpecies)-1):
        lot.append(quadOp(Symbol('t',real=True),los[tik1],v[0],s[0],True,los[tik1+1],v[0],s[0],False))
    
    if isSuper:
        for tik1 in range(0,len(listOfSpecies)-1):
            lot.append(quadOp(Symbol('Delta',real=True),los[tik1],v[0],s[0],False,los[tik1+1],v[0],s[0],False))

    for tik in range(1,len(v)):
        for tik1 in range(0,len(listOfSpecies)):
            for tik2 in range(0,len(listOfSpecies)):
                if True:#tik2 > tik1 or tik1 == tik2:
                    lot.append(quadOp(Symbol('t',real=True),los[tik1],v[0],s[0],True,los[tik2],v[tik],s[tik],False))        
                    if isSuper:
                        lot.append(quadOp(Symbol('Delta',real=True),los[tik1],v[0],s[0],False,los[tik2],v[tik],s[tik],False))

    return lot

def listStaggeredKitaev(vectorSet,vectorLabels,listOfSpecies,isSuper):

    v = vectorSet
    s = vectorLabels
    los = listOfSpecies
    lot = []
        
    for tik in range(0,len(listOfSpecies)):
	lot.append(quadOp(Symbol('mu',real=True),los[tik],v[0],s[0],True,los[tik],v[0],s[0],False))

    lot.append(quadOp(Symbol('t0',real=True),los[0],v[0],s[0],True,los[1],v[0],s[0],False))
    lot.append(quadOp(Symbol('t1',real=True),los[1],v[0],s[0],True,los[0],v[1],s[1],False))
    lot.append(quadOp(Symbol('Delta0',real=True),los[0],v[0],s[0],False,los[1],v[0],s[0],False))
    lot.append(quadOp(Symbol('Delta1',real=True),los[1],v[0],s[0],False,los[0],v[1],s[1],False))

    return lot

def listGenerate1(vectorSet, vectorLabels, listOfSpecies, isSuper):
	
    v = vectorSet
    s = vectorLabels
    los = listOfSpecies

    lot = []
    loc = []
    
    tSub = 0
    muSub = 0
    deltaSub = 0

    for tik in range(0,len(los)):
        lot.append(quadOp(Symbol('mu%d'%(muSub),real=True),los[tik],v[0],s[0],True,los[tik],v[0],s[0],False))
        loc.append(Symbol('mu%d'%(muSub),real=True))
        muSub = muSub + 1

    for tik1 in range(0,len(los)-1):
        lot.append(quadOp(Symbol('t%d'%(tSub)),los[tik1],v[0],s[0],True,los[tik1+1],v[0],s[0],False))
        loc.append(	Symbol('t%d'%(tSub)))
        tSub = tSub + 1
            

    if isSuper:
        for tik1 in range(0,len(los)-1):
            lot.append(quadOp(Symbol('Delta%d'%(deltaSub)),los[tik1],v[0],s[0],False,los[tik1+1],v[0],s[0],False))
            loc.append(Symbol('Delta%d'%(deltaSub)))
            deltaSub = deltaSub + 1

    for tik in range(1,len(v)):
        for tik1 in range(0,len(los)):
            for tik2 in range(0,len(los)):
                if True:#not tik1 == tik2:
                    lot.append(quadOp(Symbol('t%d'%(tSub)),los[tik1],v[0],s[0],True,los[tik2],v[tik],s[tik],False))        
                    loc.append(Symbol('t%d'%(tSub)))
                    tSub = tSub + 1 	

                    if isSuper:
                        lot.append(quadOp(Symbol('Delta%d'%(deltaSub)),los[tik1],v[0],s[0],False,los[tik2],v[tik],s[tik],False))
                        loc.append(Symbol('Delta%d'%(deltaSub)))
                        deltaSub = deltaSub + 1
 
    return lot

def extractCoeff(listofterms):
	
    listOfCoeff = []
    
    for tik in range(0,len(listofterms)):
        listOfCoeff.append(listofterms[tik].coeff)
    
    return listOfCoeff


















