from astropy.io import ascii
from astropy.table import Table
import numpy as np
import scipy as sp
import scipy.constants as sc
import scipy.optimize as op
import re
import sys
import optparse
import os.path
from glob import glob
import re
sys.path.append("../simulator/")
import getpriority
import Generate_Velocities
from consts import JD0, ECC, OMEGA
import SystPy

def fit_dates_vels(dates,vels,errs):
    odates = dates - np.min(dates)
    try:
        pcoeff, cov = np.polyfit(odates,vels,1,w=1.0/errs**2,cov=True)
#    print pcoeff
#    txt = "%.4f %.4f %.4f %.4f %.4f" % (pcoeff[0],pcoeff[1],np.max(odates),np.sqrt(cov[0][0]),np.sqrt(cov[1][1]))
        txt = "%.4f %.4f %.4f" % (pcoeff[0],pcoeff[1],np.max(odates))
    except:
        txt = "Cannot perform polyfit"
    return txt

def readin_velsfile(fn):
    indata = ascii.read(fn)
    names=["jd","velocity","int unc","I2 counts","phases"]
    for i in range(0,len(names)):
        oldn = "col%d" % (i+1)
        newn = names[i]
        indata.rename_column(oldn,newn)
    good = (indata['int unc'] > 0) & (indata['int unc'] < 100)
    return indata[good]

def parse_options():

    parser = optparse.OptionParser()
    parser.add_option("-a","--all",dest="all",default=False,action="store_true")
    parser.add_option("-i","--infile",dest="infile",default="../Datafiles/newgoogledex_sinnoise.csv")
    parser.add_option("-v","--veldir",dest="veldir",default="../VelsFiles/")
    parser.add_option("-o","--outdir",dest="outdir",default="../PlanetFitting/")        
    parser.add_option("-t","--tessapf",dest="mfn",default='../Datafiles/TESSAPF_AWMasses_prec_phase.csv')
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
            

    return snames,options.infile,options.mfn,veldir,options.outdir


def mcmc_planets(k,sname,TESSAPFdata,indices,veloff,writefit=False,veldir="../VelsFiles",outdir="../PlanetFitting"):


    initphases = TESSAPFdata['phase'][indices]
    vmag = TESSAPFdata['vmag'][indices]
    rplanets = TESSAPFdata['rplanet'][indices]
    mplanets = TESSAPFdata['true_mass'][indices]
    planetids = TESSAPFdata['Index'][indices]

    fn = sname + ".vels"
    indata = readin_velsfile(os.path.join(veldir,fn))
    
    m = SystPy.matrix_to_array(SystPy.getResidualMatrix(k))
    txt = fit_dates_vels(m[:,0],m[:,1],m[:,2])
    otxt = "%s %f" % (txt, k.getRMS())

    if writefit:
        outname = "%s_mcmc.fit" % (sname)
        outname = os.path.join(outdir,outname)
        outfp = open(outname,"w")
        outfp.write(otxt+"\n")


    kl = SystPy.MCMC(k,nChains=4, rStop=1.05)
    bsname = os.path.join(outdir,sname + ".kl")
    SystPy.KLSave(kl, bsname)
    stds = kl.getElementsStats(SystPy.K_STAT_STDDEV)
    meds = kl.getElementsStats(SystPy.K_STAT_MEDIAN)
    mads = 1.4826*kl.getElementsStats(SystPy.K_STAT_MAD)    
    for i in range(0,len(initphases)):
        fperiod = k.getElement(i+1,SystPy.K_PER)
        phases = getpriority.compute_currentphase(indata['jd'],fperiod,initphases[i])

        K = meds[i+1,SystPy.K_SEMIAMP]
        err_K = mads[i+1,SystPy.K_SEMIAMP]

        binedges = [getpriority.EDGE1,getpriority.EDGE2,getpriority.EDGE3,getpriority.EDGE4,getpriority.EDGE5]
        phasebins = np.digitize(phases,binedges,right=True)
        ninbins = []
        for nbin in range(0,len(binedges)):
            cpb = phasebins[phasebins == nbin]
            ninbins.append(len(cpb))

        if writefit:
            ostr = "%f %f %f %.4g %.4g %.4g %d" % (fperiod,K,err_K,k.getElement(i+1,SystPy.K_MASS)*317.8281,mads[i+1,SystPy.K_MASS]*317.8281,rplanets[i],len(phases))
            # period (days) K (m/s) err_K (m/s) planet mass (M_earth) error (M_earth) R (R_earth) #vs #phase1 #phase2 #phase3 #phase4 #phase5
            bstr = ""
            for i in range(0,len(ninbins)):
                bstr += " %d" % ninbins[i]
            bstr += "\n"
            outfp.write(ostr)
            outfp.write(bstr)

#    kernelname = "%s" % (sname)
#    kernelname = os.path.join(outdir,kernelname)
#    k.save(kernelname)
    if writefit:
        outfp.close()
    return



snames,gdfn,mfn,veldir,outdir = parse_options()

gd = ascii.read(gdfn)
TESSAPFdata=ascii.read(mfn,format='csv')
for sname in snames:

    sfn = sname + ".sys"
    vfn = sname + ".vels"

    kernelname = "%s" % (sname)
    kernelname = os.path.join(outdir,kernelname)
    if os.path.exists(kernelname):
        k = SystPy.loadKernel(kernelname)
    else:
        continue

    planets, = np.where((TESSAPFdata['star_names'] == sname) & (TESSAPFdata['detected'] == "TRUE"))
        # planet masses are in earth masses, Systemic likes Jupiters
        #mearth = 5.9722e24
        #mjup = 1.898e27
    mratio = 317.8281
    TESSAPFdata['est_mass'] /= mratio

    print "%s: RMS of fit %f" % (sname,k.getRms())

    mcmc_planets(k,sname,TESSAPFdata,planets,gd['vel_offset'][gd['starname'] == sname],veldir=veldir,outdir=outdir)

    bvfn = "%sbinned" % (sname)
    kernelname = os.path.join(outdir,bvfn)
    k = SystPy.loadKernel(kernelname)
    ddir =veldir

    print "%s: RMS of fit %f" % (bvfn,k.getRms())

    mcmc_planets(k,bvfn,TESSAPFdata,planets,gd['vel_offset'][gd['starname'] == sname],writefit=True,veldir=veldir,outdir=outdir)
        
