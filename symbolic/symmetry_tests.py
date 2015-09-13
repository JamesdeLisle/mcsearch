from sympy import *

def getPauli(number):
    
    if number == 0:
        return Matrix([[1,0],[0,1]])
    elif number == 1:
        return Matrix([[0,1],[1,0]])
    elif number == 2:
        return Matrix([[0,-1j],[1j,0]])
    elif number == 3:
        return Matrix([[1,0],[0,-1]])

def testPH(CP,H):

    return CP * conjugate(H) * CP.adjoint() + H.subs([(Symbol('p0',real=True),-Symbol('p0',real=True)),(Symbol('p1',real=True),-Symbol('p1',real=True)),(Symbol('p2',real=True),-Symbol('p2',real=True))])

def testTR(CT,H):
    
    return CT * conjugate(H) * CT.adjoint() - H.subs([(Symbol('p0',real=True),-Symbol('p0',real=True)),(Symbol('p1',real=True),-Symbol('p1',real=True)),(Symbol('p2',real=True),-Symbol('p2',real=True))])

def constructPHInvH(CP,H):

    return H.subs([(Symbol('p0',real=True),-Symbol('p0',real=True)),(Symbol('p1',real=True),-Symbol('p1',real=True)),(Symbol('p2',real=True),-Symbol('p2',real=True))]) - CP * conjugate(H) * CP.adjoint()

def constructTRInvH(TR,H):

    return H.subs([(Symbol('p0',real=True),-Symbol('p0',real=True)),(Symbol('p1',real=True),-Symbol('p1',real=True)),(Symbol('p2',real=True),-Symbol('p2',real=True))]) + TR * conjugate(H) * TR.adjoint()
