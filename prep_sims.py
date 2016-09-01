import os
import shutil
from glob import glob

fns = glob("./201?-??-??.simout")
for fn in fns:
    try:
        os.unlink(fn)
    except:
        print "File %s no longer exists or you lack the permissions" % (fn)
        

otherfns = ["./newgoogledex.csv"]
for fn in fns:
    try:
        os.unlink(fn)
    except:
        print "File %s no longer exists or you lack the permissions" % (fn)

shutil.copyfile("../Datafiles/newgoogledex_sinnoise.csv","./newgoogledex.csv")
        
fns = glob("../SimFiles/*.sim")
for fn in fns:
    try:
        os.unlink(fn)
    except:
        print "File %s no longer exists or you lack the permissions" % (fn)
        

fns = glob("../VelsFiles/*.vels")
for fn in fns:
    try:
        os.unlink(fn)
    except:
        print "File %s no longer exists or you lack the permissions" % (fn)
