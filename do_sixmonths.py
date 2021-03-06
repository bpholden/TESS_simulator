import subprocess
import argparse
import shutil
import shlex
import os
from glob import glob
parser = argparse.ArgumentParser()
parser.add_argument("scheme")
parser.add_argument("old_googledex")
parser.add_argument("prefix")
parser.add_argument("startdate")
parser.add_argument("enddate")
parser.add_argument("-s","--seed",help="optional seed",type = int)
args = parser.parse_args()

# need old_googledex new_googledex startdate enddate

scheme = args.scheme
old_googledex = args.old_googledex
prefix  = args.prefix
startdate = args.startdate
enddate = args.enddate
if args.seed:
    seed = args.seed
else:
    seed = None

outdir = prefix + "_sinnoise_" + scheme + "_noa_singlestandard_twothirds"

new_googledex = "newgoogledex_" + prefix + ".csv"

exstr = "python sim_nights.py -p %s -i %s %s %s -d" % (scheme,old_googledex,startdate,enddate)
if seed:
    exstr += " -s %d" %(seed)
print exstr
out = subprocess.check_output(exstr,shell=True)
print out
#exstr = "python make_sim_files.py"
#subprocess.check_output(exstr,shell=True)
exstr =  "python makevels.py"
subprocess.check_output(exstr,shell=True)

#../SystPy/
exstr = "python fit_TESS_APF.py -i %s -a" % (os.path.join("../simulator/",old_googledex))
out = subprocess.check_output(exstr,shell=True,cwd="../SystPy/")
print out
# ../simulator/
if os.path.exists(new_googledex):
    prev_googledex = new_googledex + ".prev"
    os.rename(new_googledex,prev_googledex)

exstr = "python TESSAPF_assess.py -i %s -o %s -p ../PlanetFitting/" %(old_googledex,new_googledex)
out = subprocess.check_output(exstr,shell=True)
fn_assess = "TESSAPF_assess_" + prefix + ".out"
fp = open(fn_assess,"w")
fp.write(out)
fp.close()
print out
exstr  = "python fillin_phasebins.py -i %s " % (new_googledex)
try:
    output = subprocess.check_output(exstr,shell=True)
    print output
except:
    pass 
for id in ["../SimFiles","../VelsFiles/","../PlanetFitting/"]:
    foutdir = os.path.join(id,outdir)
    if not os.path.exists(foutdir):
        os.mkdir(foutdir)
    else:
        print "%s already exists!" % (foutdir)

vels = glob("../VelsFiles/TESS*")
for vel in vels:
    shutil.move(os.path.join("../VelsFiles/",vel),os.path.join("../VelsFiles/",outdir))

pfs = glob("../PlanetFitting/TESS*")
for pf in pfs:
    shutil.move(os.path.join("../PlanetFitting/",pf),os.path.join("../PlanetFitting/",outdir))

shutil.move(fn_assess,os.path.join("../PlanetFitting/",outdir))
shutil.move(old_googledex,os.path.join("../SimFiles/",outdir))

sims = glob("../SimFiles/*.sim")
for sim in sims:
    shutil.copy(os.path.join("../SimFiles/",sim),os.path.join("../SimFiles/",outdir))
    
