from sympy import  *
from IPython.display import display
import copy

class quadOp:           

    #################################################
    # Set symbol based on current internal variables 
    def setSymbol(self):
        dag = "\dagger"
        
        if self.DAGboola == True:
            self.symbol1 = Symbol("{{%s}_{%s}^{%s}}" % (self.typea,self.indexaSym,dag), commutative=False)
        else:
            self.symbol1 = Symbol("{{%s}_{%s}}" % (self.typea,self.indexaSym), commutative=False)
        
        if self.DAGboolb == True:
            self.symbol2 = Symbol("{{%s}_{%s}^{%s}}" % (self.typeb,self.indexbSym,dag), commutative=False)
        else:
            self.symbol2 = Symbol("{{%s}_{%s}}" % (self.typeb,self.indexbSym), commutative=False)
        
        self.symbol = self.coeff*self.symbol1*self.symbol2
    #################################################

    #################################################
    # Initate internal variables of term
    def __init__(self, coeff, typea, indexa, indexaSym, DAGboola, typeb, indexb, indexbSym, DAGboolb):        
        
        self.coeff = coeff
        self.typea = typea
        self.typeb = typeb
        self.indexa = indexa
        self.indexb = indexb
        self.indexaSym = indexaSym
        self.indexbSym = indexbSym
        self.DAGboola = DAGboola
        self.DAGboolb = DAGboolb
        self.setSymbol()
    #################################################

    #################################################
    # Perform Hermitian conjugate of term
    def HC(self):
        
        if self.DAGboola == self.DAGboolb:
            self.DAGboola = not self.DAGboola
            self.DAGboolb = not self.DAGboolb
        
        temptypea = self.typea
        self.typea = self.typeb
        self.typeb = temptypea
        tempindexa = self.indexa
        self.indexa = self.indexb
        self.indexb = tempindexa
        tempindexaSym = self.indexaSym
        self.indexaSym = self.indexbSym
        self.indexbSym = tempindexaSym
        self.coeff = conjugate(self.coeff)    
        self.setSymbol()
    #################################################

    #################################################
    # Perform Fourier transform of term    
    def FT(self):        

        p = symbols('p0:%d' % (len(self.indexa)), real=True)
        q = Symbol('p',real=True)
        k = 0
        
        if self.indexb == zeros(len(self.indexa),1):
            for tik in range(0,len(self.indexa)):
                k = k + copy.deepcopy(p[tik])*copy.deepcopy(self.indexa[tik])
                
        else:
             for tik in range(0,len(self.indexb)):
                k = k + copy.deepcopy(p[tik])*copy.deepcopy(self.indexb[tik])    
                    
        if self.indexa-self.indexb == zeros(len(self.indexa),1):
            pp = 0;
        else:
            if self.indexa != zeros(len(self.indexa),1):
                pp = 1
            else:
                pp = -1 
        
        

        if not self.DAGboola and not self.DAGboolb:   
            self.indexbSym = -q
            self.coeff = self.coeff * exp(I*pp*k)
        
        if self.DAGboola and self.DAGboolb:
            self.indexbSym = -q
            self.coeff = self.coeff * exp(-I*pp*k)
        
        if self.DAGboola and not self.DAGboolb:
            self.indexbSym = q
            self.coeff = self.coeff * exp(-I*pp*k)
        
        if not self.DAGboola and self.DAGboolb:
            self.indexbSym = q
            self.coeff = self.coeff * exp(I*pp*k)
        
        self.indexaSym = q
        self.setSymbol()
    ################################################# 
    
    #################################################
    # Takes term from p to -p
    def minusP(self):
       
        q = Symbol('p',real=True)
        qq = -Symbol('p',real=True)        

        if self.indexaSym == Symbol('p',real=True):
            self.indexaSym = -Symbol('p',real=True)
        else:
            self.indexaSym = Symbol('p',real=True)
        if self.indexbSym == Symbol('p',real=True):
            self.indexbSym = -Symbol('p',real=True)
        else:
            self.indexbSym = Symbol('p',real=True)
        
        self.coeff = simplify(self.coeff.subs([(Symbol('p0',real=True),-Symbol('p0',real=True)),(Symbol('p1',real=True),-Symbol('p1',real=True)),(Symbol('p2',real=True),-Symbol('p2',real=True))]))
        self.coeff = simplify(self.coeff) 
        self.setSymbol()
    #################################################

    def commute(self):
        pass        



