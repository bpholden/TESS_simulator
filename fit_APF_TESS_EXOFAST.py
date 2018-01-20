from astropy.io import ascii
from astropy.table import Table
import numpy as np
import scipy as sp
import pandas as pd
import scipy.constants as sc
import scipy.optimize as op
import matplotlib.pyplot as plt
import re
import sys
import optparse
import copy
import os.path
from glob import glob
import re
sys.path.append("../simulator/")
import getpriority
import Generate_Velocities
from consts import JD0, ECC, OMEGA
import ConfigParser


def parse_options():

    parser = optparse.OptionParser()
    parser.add_option("-a","--all",dest="all",default=False,action="store_true")
    parser.add_option("-i","--infile",dest="infile",default="../Datafiles/newgoogledex.csv")
    parser.add_option("-v","--veldir",dest="veldir",default="../VelsFiles/")
    parser.add_option("-o","--outdir",dest="outdir",default="../PlanetFitting/")        
    parser.add_option("-t","--tessapf",dest="mfn",default='../Datafiles/TESSAPF_AWMasses_prec_phase.csv')
    parser.add_option("-e","--extra",dest="extra",default=False)
    (options, args) = parser.parse_args()
    if len(args) < 1 and not options.all:
        print "needs a star name"
        sys.exit()
    snames = []

    veldir = options.veldir + "/"
    
    if options.all:
        files = glob(veldir+"TESSAPF*.vels")
        snames = []
        for f in files:
            m = re.search("(TESSAPF\d+)\.vels",f)
            if m:
                snames.append(m.group(1))
    else:
        for a in args:
            if re.match("\A\d+\Z",a):
                sname = "TESSAPF%s" % (a)
            else:
                sname = a
            vname = os.path.join(veldir,sname + ".vels")
            if not os.path.exists(vname):
                print "%s does not exist" %(sname)
                sys.exit()
            snames.append(sname)

    if not os.path.isdir(options.outdir):
        try:
            os.mkdir(options.outdir)
        except Exception as e:
            print "cannot make dir %s: %s" % (options.outdir,e)
            sys.exit()
            

    return snames,options.infile,options.mfn,veldir,options.outdir,options.extra

def initialize_model(sname,planets,TESSAPFdata,addextra=False):
    time_base = 2458362.5

    priorname = sname + ".priors"
    priorfp = open(priorname,"w+")
    fitname = "fit_" + sname + ".pro"
    fitfp = open(fitname,"w+")

    fitfp.write("pro fit_" + sname + ",debug=debug\n")
    fitfp.write("\nmaxsteps=100000\nnthin=2\n")

    nplan = len(planets)
    if addextra:
        nplan = nplan + 1

    rvpath = "%s.vels" % (sname)
    execstr = "exofastv2, nplanets=%d,rvpath='%s',priorfile='%s',prefix='%s'," % (nplan,rvpath,priorname,sname)
    execstr = execstr + "maxsteps=maxsteps,nthin=nthin,"
    trues = ""
    falses = ""
    for i in range(0,nplan):
        trues += "1,"
        falses +=  "0,"
    trues = trues.rstrip(",")
    falses = falses.rstrip(",")
    execstr = execstr + "circular=[%s],fitrv=[%s]" % (trues,trues)
    execstr = execstr + ",debug=debug"
    fitfp.write(execstr + "\n")
    fitfp.write("\nend\n")
    fitfp.close()
    n = 0
    ks = Generate_Velocities.calc_K(TESSAPFdata,sname)
    
    priorfp.write("mstar %f 0\n" % (TESSAPFdata['mstar'][planets[0]]))  
    for idx in planets:
        c_period = TESSAPFdata['period'][idx]
        priorfp.write("period_%d %f 0\n" % (n,c_period))
        priorfp.write("k_%d %f 15.\n" % (n,float(ks[n])))
        tpstart = time_base + c_period * (TESSAPFdata['phase'][idx] + 0.25)
        priorfp.write("tc_%d %f 15.\n" % (n,tpstart))
#        priorfp.write("cosi_%d 0. 0." % (idx))                         
        n = n + 1

    if addextra:

        priorfp.write("period_%d 30 15\n" % (n))
        priorfp.write("k_%d 2.\n" % (n))
        tpstart = time_base + 15.
        priorfp.write("tc_%d %f\n" % (n,tpstart))
        priorfp.write("cosi_%d 0. 0.\n" % (idx))                         
        
    priorfp.close()

    execstr = "idl -e '%s'" %(fitname)

    return execstr


if __name__ == "__main__":

    snames,gdfn,mfn,veldir,outdir,addextra = parse_options()

    gd = ascii.read(gdfn)
    TESSAPFdata=ascii.read(mfn,format='csv')
    # planet masses are in earth masses, Systemic likes Jupiters
    #mearth = 5.9722e24
    #mjup = 1.898e27
    mratio = 317.83
    TESSAPFdata['est_mass'] /= mratio
            
    for sname in snames:

        planets, = np.where((TESSAPFdata['star_names'] == sname) & (TESSAPFdata['detected'] == "TRUE"))
        estr = initialize_model(sname,planets,TESSAPFdata,addextra=addextra)
        print ( estr )
