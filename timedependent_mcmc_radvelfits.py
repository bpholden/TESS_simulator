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
import radvel
import ConfigParser
import Radvel_MCMC

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


if __name__ == "__main__":

    snames,gdfn,mfn,veldir,outdir,addextra = parse_options()

    gd = ascii.read(gdfn)
    TESSAPFdata=ascii.read(mfn,format='csv')
    # planet masses are in earth masses, Systemic likes Jupiters
    #mearth = 5.9722e24
    #mjup = 1.898e27
    mratio = 317.8281
    TESSAPFdata['est_mass'] /= mratio
            
    for sname in snames:

        sfn = sname + ".sys"
        vfn = sname + ".vels"
        newvels = True

        invels = Radvel_MCMC.readin_velsfile(os.path.join(veldir,vfn))
        ddir =veldir

        current_nvels =0
        for i in range(0,12):

            deltas = (invels['time'] - invels['time'].min() )< 90*(i+1)
            if len(invels['time'][deltas]) <= current_nvels + 5:
                continue

            current_nvels = len(invels['time'][deltas])
 
            planets, = np.where((TESSAPFdata['star_names'] == sname) & (TESSAPFdata['detected'] == "TRUE"))

            mod = Radvel_MCMC.initialize_model(planets,TESSAPFdata,addextra=addextra)
            like = Radvel_MCMC.make_like(mod,invels[deltas],planets,addextra=addextra)
            post = Radvel_MCMC.init_posterior(like,planets,addextra=addextra)
            post = radvel.fitting.maxlike_fitting(post, verbose=False)

            chains,Ks,err_Ks,Ms,err_Ms = Radvel_MCMC.mcmc_planets(post,outdir,sname,TESSAPFdata['mstar'][planets],addextra=addextra)
            Radvel_MCMC.writefit(sname,outdir,TESSAPFdata['period'][planets],Ks,err_Ks,Ms,err_Ms,TESSAPFdata['rplanet'][planets],len(invels))
