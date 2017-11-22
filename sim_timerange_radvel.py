import subprocess
import argparse
import shutil
import shlex
import sys
import os
from glob import glob


allowed_phases = ["nights","vels","fit","cat"]
phasehelp = "optional phase, allowed: " + " ".join(allowed_phases)
parser = argparse.ArgumentParser()
parser.add_argument("scheme")
parser.add_argument("prefix")
parser.add_argument("startdate")
parser.add_argument("enddate")
parser.add_argument("-s","--seed",help="optional seed",type = int)
parser.add_argument("-p","--phase",help=phasehelp)
args = parser.parse_args()

scheme = args.scheme
prefix  = args.prefix
startdate = args.startdate
enddate = args.enddate
if args.seed:
    seed = args.seed
else:
    seed = None

if args.phase:
    phase = args.phase
else:
    phase = None


if phase and phase not in allowed_phases:
    print "%s is not an allowed phase" % (phase)
    sys.exit()
    
outdir = prefix + "_" + scheme + "_noa_twothirds"
pfdir  = os.path.join("..","PlanetFitting",outdir)
simdir = os.path.join("..","SimFiles",outdir)
veldir = os.path.join("..","VelsFiles",outdir)

old_googledex = "newgoogledex.csv"
fpold_googledex = os.path.join(outdir,old_googledex)

new_googledex = "newgoogledex_" + prefix + ".csv"
new_googledex = os.path.join(outdir,new_googledex)
    
if phase == "nights" or phase == None:
    if not os.path.isdir(outdir):
        try:
            os.mkdir(outdir)
        except Exception as e:
            print "cannot make %s: %s" % (outdir,e)

    shutil.copyfile("../Datafiles/newgoogledex_sinnoise.csv",fpold_googledex)


    exstr = "python sim_nights.py -p %s -i %s %s %s -d -o %s" % (scheme,old_googledex,startdate,enddate,outdir)
    if seed:
        exstr += " -s %d" %(seed)
    print exstr
    out = subprocess.check_output(exstr,shell=True)
    print out
#exstr = "python make_sim_files.py"
#subprocess.check_output(exstr,shell=True)


if phase == "vels" or phase == None:

    exstr =  "python makevels.py -i %s -o %s " % (simdir,veldir)
    subprocess.check_output(exstr,shell=True)

    
if phase == "fit" or phase == None:
    #../SystPy/
    exstr = "python fit_APF_TESS_radvel.py -i %s -v %s -o %s -a" % (os.path.join("../simulator/",fpold_googledex),veldir,pfdir)
    out = subprocess.check_output(exstr,shell=True,cwd="../SystPy/")
    print out
    # ../simulator/


if phase == "cat" or phase == None:

    outcat = "mass_radius" + prefix + "_" + scheme + ".out"
    
    exstr = "python make_radvel_catalogs.py -i %s -o %s" %(pfdir,outcat)
    out = subprocess.check_output(exstr,shell=True)
    print out
