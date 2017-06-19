import sqlite3 as lite
import sys
import os
import re

import Qdata as QD

#Grep for a couple things to determine whether this is actually a Qchem output file
#This does not assume that the calculation was actually successful
def checkQcOut(f):
    counter = 0
    str1 = 'Welcome to Q-Chem'
    str2 = 'A Quantum Leap Into The Future Of Chemistry'
    with open(f) as qcout:
        head = [next(qcout) for line in range(10)]
    checkstr1 = [str1 in line for line in head]
    checkstr2 = [str2 in line for line in head]
    if True in checkstr1 and True in checkstr2: return True
    else: return False

            


#convert a QchemData class object into data in the sqlite database
def sqliteQcOut(qcout,database):
    pass



con = lite.connect(sys.argv[1]) #connect with database
datadir = sys.argv[2]
dirpath = os.getcwd()

for filename in os.listdir(datadir):
    fpath = '/'.join([dirpath,datadir,filename])
    if checkQcOut(fpath):
        qcout = QD.Qdata(fpath)
        sqliteQcOut(qcout,con)
        







