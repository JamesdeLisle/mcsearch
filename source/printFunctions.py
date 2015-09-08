from __future__ import print_function
import random
from numpy import real, dot, absolute, array, conjugate, linalg
import matplotlib.pyplot as plt
import time
from tabulate import tabulate
import curses

def initialiseDisplay():

    stdscr = curses.initscr()
    curses.noecho()
    curses.cbreak()
    
    return stdscr

def killDisplay():

    curses.nocbreak()
    curses.echo()
    curses.endwin()

def printTable(stdscr, step,step_maximum,temperature,costFunction_current,costFunction_minimum):
    
    out = [['Step','%d/%d' %(step,step_maximum)],['Inverse Temperature','%f' % (temperature)],['Cost Function','%f' % (costFunction_current)],['Cost Function Minimum','%f' % (costFunction_minimum)]]
    stdscr.addstr(30,0,tabulate(out,tablefmt='grid'))
    stdscr.refresh()

def printLambdaFlags(stdscr,activeFlags):
    
    out = []
    variables = ['mu0','mu1','t0','Delta0','t1','Delta1','t2','Delta2','t3','Delta3','t4','Delta4']

    for tik in range(12):
        out.append([variables[tik],activeFlags[0][tik],activeFlags[1][tik]])
        stdscr.addstr(0,0,tabulate(out,headers=['Variable','Real','Imaginary'],tablefmt='grid'))
        stdscr.refresh()

def printLambda(stdscr,Lambda):
    
    out = []
    variables = ['mu0','mu1','t0','Delta0','t1','Delta1','t2','Delta2','t3','Delta3','t4','Delta4']

    for tik in range(12):
        out.append([variables[tik],round(Lambda[0][tik],3),round(Lambda[1][tik],3)])
        stdscr.addstr(0,0,tabulate(out,headers=['Variable','Real','Imaginary'],tablefmt='grid',stralign="center"))
        stdscr.refresh()

def printLambdaFinal(Lambda_current):

    stdscr = initialiseDisplay()
    printLambdaFlags(stdscr,Lambda_current)
    y,x = curses.getsyx()
    stdscr.addstr(y+2,0,"Press 'q' to return to the console")
    stdscr.refresh() 
    key = ''
    while key != ord('q'):
        key = stdscr.getch()
    stdscr.clear()
    killDisplay()

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

def printRunToFile(costFunc_min,Lambda_min,Lambda,costFunc,CTR,temperature_maximum,momentum_discretisation):
    

    f = open(time.strftime('%d-%m-%y')+' '+time.strftime('%H:%M')+' '+'Max: %d, Mom: %d, Cost_Min: %f.dat'%(temperature_maximum,momentum_discretisation,costFunc_min),'a')
    TextList = ['m0','m1','t0','D0','t1','D1','t2','D2','t3','D3','t4','D4']
    print('%-------------------------------------------------%',file=f)
    print(time.strftime('%H:%M:%S'),file=f)
    print(time.strftime('%d\%m\%Y'),file=f)
    print('Max Temperature = %f' % (temperature_maximum),file=f) 
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


