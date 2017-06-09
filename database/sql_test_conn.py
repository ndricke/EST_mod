#!/usr/bin/python
# -*- coding: utf-8 -*-
 
import sqlite3 as lite
import sys
 
 
con = lite.connect('system.db')
 
with con:    
 
    cur = con.cursor()    
    cur.execute("SELECT users.name, jobs.profession FROM jobs INNER JOIN users ON users.ID = jobs.uid")
 
    rows = cur.fetchall()
 
    for row in rows:
        print row
