from dmetparse import DmetParse
#import Pandas as pd
import os
import sys

head = "lattype,lat,part,m,E,E_cen\n"

indir = sys.argv[1]
path = indir + "/"
print path
outfile = sys.argv[2]

try:
    if os.stat(outfile).st_size > 0: pass
    else:
        with open(outfile,'w') as f:
            f.write(head)
except OSError:
        with open(outfile,'w') as f:
            f.write(head)

for filename in os.listdir(indir):
    fp = path+filename
    if os.path.isfile(fp) and "DMRG" not in filename: 
        print "Parsing: " + filename
        dm = DmetParse(fp)
        if dm.E == None or dm.E_cen == None: print fp + " did not complete successfully!"
        else: dm.dataWrite(outfile)

    else:
        print "Was reading either a directory or DMRG job; continuing..."
        pass




