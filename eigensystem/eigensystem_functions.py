import numpy

def getKey(item):

    return item[0]

def sortEigensystem(eigenvalues,eigenvectors):
        
    sort = sorted(zip(eigenvalues,eigenvectors),key=getKey)
    out = zip(*sort)
    evals_sorted = numpy.array(out[0])
    evec_sorted = numpy.array(out[1])

    return evals_sorted,evec_sorted  

