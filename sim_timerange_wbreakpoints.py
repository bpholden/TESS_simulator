import subprocess
import argparse
import shutil
import shlex
import sys
import os
from glob import glob
from datetime import datetime,timedelta

def parsedate(datestr):
    cdate = datetime.strptime(datestr,"%Y/%m/%d")
    return cdate

def strdate(dtobj):
    strdate = "%d/%d/%d" % (dtobj.year,dtobj.month,dtobj.day)
    return strdate

parser = argparse.ArgumentParser()
parser.add_argument("scheme")
parser.add_argument("prefix")
parser.add_argument("startdate")
parser.add_argument("enddate")
parser.add_argument("-s","--seed",help="optional seed",type = int)
args = parser.parse_args()

scheme = args.scheme
prefix  = args.prefix
startdate = args.startdate
enddate = args.enddate
if args.seed:
    seed = args.seed
else:
    seed = None

    
outdir = prefix + "_" + scheme + "_wa_twothirds"
pfdir  = os.path.join("..","PlanetFitting",outdir)
simdir = os.path.join("..","SimFiles",outdir)
veldir = os.path.join("..","VelsFiles",outdir)

old_googledex = "newgoogledex.csv"
fpold_googledex = os.path.join(outdir,old_googledex)

new_googledex = "newgoogledex_" + prefix + ".csv"
new_googledex = os.path.join(outdir,new_googledex)

currentdt = parsedate(startdate)
enddt = parsedate(enddate)
zerotd = timedelta()

while enddt - currentdt  > zerotd:
    

    if not os.path.isdir(outdir):
        try:
            os.mkdir(outdir)
        except Exception as e:
            print "cannot make %s: %s" % (outdir,e)

    shutil.copyfile("../Datafiles/newgoogledex_sinnoise.csv",fpold_googledex)

    cstartdate = strdate(currentdt)
    newmonth = currentdt.month+6
    newyear = currentdt.year
    if newmonth > 12:
        newmonth %= 12
        newyear += 1
    currentdt = datetime(newyear,newmonth,currentdt.day)
    cenddate = strdate(currentdt)
        
    exstr = "python sim_nights.py -p %s -i %s %s %s -d -o %s" % (scheme,old_googledex,cstartdate,cenddate,outdir)
    if seed:
            exstr += " -s %d" %(seed)
            print exstr
            
    out = subprocess.check_output(exstr,shell=True)
    print out
    #exstr = "python make_sim_files.py"
    #subprocess.check_output(exstr,shell=True)



    exstr =  "python makevels.py -i %s -o %s " % (simdir,veldir)
    subprocess.check_output(exstr,shell=True)

    

    #../SystPy/
    exstr = "python fit_TESS_APF.py -i %s -v %s -o %s -a" % (os.path.join("../simulator/",fpold_googledex),veldir,pfdir)
    out = subprocess.check_output(exstr,shell=True,cwd="../SystPy/")
    print out
    # ../simulator/


    if os.path.exists(new_googledex):
        prev_googledex = new_googledex + ".prev"
        os.rename(new_googledex,prev_googledex)

    exstr = "python TESSAPF_assess.py -i %s -o %s -p %s" %(fpold_googledex,new_googledex,pfdir)
    out = subprocess.check_output(exstr,shell=True)
    print out
    fn_assess = "TESSAPF_assess_" + prefix + ".out"
    fn_assess = os.path.join(pfdir,fn_assess)
    fp = open(fn_assess,"w")
    fp.write(out)
    fp.close()

    exstr  = "python fillin_phasebins.py -i %s " % (new_googledex)
    try:
        output = subprocess.check_output(exstr,shell=True)
        print output
    except:
        pass 

    
outcat = "mass_radius" + prefix + "_" + scheme + ".out"
    
exstr = "python count_planets.py -i %s -o %s" %(pfdir,outcat)
out = subprocess.check_output(exstr,shell=True)
print out
