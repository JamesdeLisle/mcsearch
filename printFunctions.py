from __future__ import print_function
import random
from numpy import real, dot, absolute, array, conjugate, linalg
from math import cos,sin,exp
import cmath as cm
import matplotlib.pyplot as plt
import time


def printLambdaToConsole(Lambda):
    
    TextList = ['mu0','mu1','t0','Delta0','t1','Delta1','t2','Delta2','t3','Delta3','t4','Delta4']
    for	tik in range(0,12):
        if Lambda[0][tik] == 1.0 and Lambda[1][tik] != 1.0:
            print('%s\t = ANY + %fi' % (TextList[tik],Lambda[1][tik]))
        elif Lambda[0][tik] != 1.0 and Lambda[1][tik] == 1.0:
            print('%s\t = %f + ANYi' % (TextList[tik],Lambda[0][tik],))
        elif Lambda[1][tik] == 1.0 and Lambda[0][tik] == 1.0:
            print('%s\t = ANY + ANYi' % (TextList[tik]))
        else:
            print('%s\t = %f + %fi' % (TextList[tik],Lambda[0][tik],Lambda[1][tik]))


def printRunToFile(costFunc_min,Lambda_min,Lambda,costFunc,CTR,temperature_maximum,temperature_increment,momentum_discretisation):
    

    f = open(time.strftime('%d-%m-%y')+' '+time.strftime('%H:%M')+' '+'Max: %d, Step: %f, Mom: %d, Cost_Min: %f.dat'%(temperature_maximum,temperature_increment,momentum_discretisation,costFunc_min),'a')
    TextList = ['m0','m1','t0','D0','t1','D1','t2','D2','t3','D3','t4','D4']
    print('%-------------------------------------------------%',file=f)
    print(time.strftime('%H:%M:%S'),file=f)
    print(time.strftime('%d\%m\%Y'),file=f)
    print('Max Temperature = %f' % (temperature_maximum),file=f)
    print('Temperature step = %f' % (temperature_increment),file=f)
    print('Momentum discretisation = %d' % (momentum_discretisation),file=f)
    print('\n',file=f)
    print('___________________________________',file=f)
    print('Final CostFunction = %f' % (costFunc),file=f)
    print('\n',file=f)
    for	tik in range(0,12):
        if Lambda[0][tik] == 1.0 and Lambda[1][tik] != 1.0:
            print('%s\t = ANY + %fi' % (TextList[tik],Lambda[1][tik]),file=f)
        elif Lambda[0][tik] != 1.0 and Lambda[1][tik] == 1.0:
            print('%s\t = %f + ANYi' % (TextList[tik],Lambda[0][tik],),file=f)
        elif Lambda[1][tik] == 1.0 and Lambda[0][tik] == 1.0:
            print('%s\t = ANY + ANYi' % (TextList[tik]),file=f)
        else:
            print('%s\t = %f + %fi' % (TextList[tik],Lambda[0][tik],Lambda[1][tik]),file=f)
    print('___________________________________',file=f)
    print('\n',file=f)
    print('___________________________________',file=f)
    print('Minimum CostFunction = %f' % (costFunc_min),file=f)
    print('\n',file=f)
    for	tik in range(0,12):
        if Lambda_min[0][tik] == 1.0 and Lambda_min[1][tik] != 1.0:
            print('%s\t = ANY + %fi' % (TextList[tik],Lambda_min[1][tik]),file=f)
        elif Lambda_min[0][tik] != 1.0 and Lambda[1][tik] == 1.0:
            print('%s\t = %f + ANYi' % (TextList[tik],Lambda_min[0][tik],),file=f)
        elif Lambda_min[1][tik] == 1.0 and Lambda_min[0][tik] == 1.0:
            print('%s\t = ANY + ANYi' % (TextList[tik]),file=f)
        else:
            print('%s\t = %f + %fi' % (TextList[tik],Lambda_min[0][tik],Lambda_min[1][tik]),file=f)
    print('___________________________________',file=f)
    print('\n',file=f)
    print('CTR = ',file=f)
    for tik in range(0,4):
        print('|\t %d+%di\t %d+%di\t %d+%di\t %d+%di\t |' % (CTR[tik][0].real,CTR[tik][0].imag,CTR[tik][1].real, \
                CTR[tik][1].imag,CTR[tik][2].real,CTR[tik][2].imag,CTR[tik][3].real,CTR[tik][3].imag),file=f)
    print('\n',file=f)
    print('%-------------------------------------------------%',file=f)
    f.close()


