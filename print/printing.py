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


