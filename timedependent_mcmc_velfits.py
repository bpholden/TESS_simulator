from __future__ import print_function

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
import fit_TESS_APF 


def bin_phase_dates(dates,phases,velocities,uncs,i2counts):

    dailydates = []
    dailyvels = []
    dailyi2sum = []
    dailyerrs = []
    dailyphases = []
    done = np.zeros_like(dates,dtype=np.bool)
    while not done.all():
        notdone, = np.where(done==False)
        cdate = dates[notdone[0]]
        cphase = phases[notdone[0]]
        ddates = dates - cdate
        dphase = phases - cphase
        curday = (done == False) & (ddates < 1.) & (ddates >= 0) & (dphase < 0.1) # same day
        curphases = phases[curday]
        curdates = dates[curday]
        curvels = velocities[curday]
        curuncs = uncs[curday]
        curi2 = i2counts[curday]
        wgts = 1.0/(curuncs*curuncs)
        dailydate = np.average(curdates)
        dailyvel = np.average(curvels,weights=wgts)
        dailyerr = np.average(curuncs)/np.sqrt(len(curuncs))
        dailyphases.append(np.average(curphases,weights=wgts))
        dailyi2 = np.sum(curi2)
#        if len(uncs) > 2:
#            dailyerr /= np.sqrt(len(uncs) - 1)
        dailydates.append(dailydate)
        dailyvels.append(dailyvel)
        dailyerrs.append(dailyerr)
        dailyi2sum.append(dailyi2)
        done[curday] = True
    return dailydates,dailyphases,dailyvels,dailyerrs,dailyi2sum


def parse_options():

    parser = optparse.OptionParser()
    parser.add_option("-a","--all",dest="all",default=False,action="store_true")
    parser.add_option("-i","--infile",dest="infile",default="../Datafiles/newgoogledex_sinnoise.csv")
    parser.add_option("-v","--veldir",dest="veldir",default="../VelsFiles/")
    parser.add_option("-o","--outdir",dest="outdir",default="../PlanetFitting/")
    parser.add_option("-f","--outfile",dest="outfn",default="../PlanetFitting/sn_v_time")            
    parser.add_option("-t","--tessapf",dest="mfn",default='../Datafiles/TESSAPF_AWMasses_prec_phase.csv')
    (options, args) = parser.parse_args()
    if len(args) < 1 and not options.all:
        print ("needs a star name")
        sys.exit()
    snames = []

    veldir = options.veldir + "/"
    binned = False
    if options.all:
        files = glob(veldir+"TESSAPF*.vels")
        snames = []
        for f in files:
            m = re.search("(TESSAPF\d+)\.vels",f)
            if m:
                snames.append(m.group(1))
            else:
                m = re.search("(TESSAPF\d+)binned.vels",f)
                binned = True
                
    else:
        for a in args:
            if re.match("\A\d+\Z",a):
                sname = "TESSAPF%s" % (a)
            else:
                sname = a
            vname = os.path.join(veldir,sname + ".vels")
            if not os.path.exists(vname):
                print ("%s does not exist" %(sname))
                sys.exit()
            snames.append(sname)

    if not os.path.isdir(options.outdir):
        try:
            os.mkdir(options.outdir)
        except Exception as e:
            print ("cannot make dir %s: %s" % (options.outdir,e))
            sys.exit()
            

    try:
        outfile = open(options.outfn,"a+")
    except Exception as e:
        print ("cannot open file %s: %s" % (options.outfn,e))
        sys.exit()
            
    return snames,options.infile,options.mfn,veldir,options.outdir,outfile,binned



def fit_planets(k,sname,initphases,vmag,rplanets,mplanets,planetids,veloff,outfile,writefit=False,veldir="../VelsFiles",outdir="../PlanetFitting"):
    
    fn = sname + ".vels"
    indata = fit_TESS_APF.readin_velsfile(os.path.join(veldir,fn))

    m = SystPy.matrix_to_array(SystPy.getResidualMatrix(k))
    txt = fit_TESS_APF.fit_dates_vels(m[:,0],m[:,1],m[:,2])
    otxt = "%s %f" % (txt, k.getRMS())
    print (otxt)

    if writefit:
        outname = "%s.fit" % (sname)
        outname = os.path.join(outdir,outname)
        outfp = open(outname,"w")
        outfp.write(otxt+"\n")

    kl = SystPy.MCMC(k,nChains = 4, rStop = 1.05)

    stds = kl.getElementsStats(SystPy.K_STAT_STDDEV)
    meds = kl.getElementsStats(SystPy.K_STAT_MEDIAN)
    mads = 1.4826*kl.getElementsStats(SystPy.K_STAT_MAD)    
    for i in range(0,len(initphases)):
        fperiod = k.getElement(i+1,SystPy.K_PER)
        phases = getpriority.compute_currentphase(indata['jd'],fperiod,initphases[i])
        

        K = meds[i+1,SystPy.K_SEMIAMP]
        err_K = mads[i+1,SystPy.K_SEMIAMP]

        ostr = "%s %f %f %f %.4g %.4g %.4g %d\n" % (sname,fperiod,K,err_K,k.getElement(i+1,SystPy.K_MASS)*317.8281,mads[i+1,SystPy.K_MASS]*317.8281,rplanets[i],len(phases))
        # period (days) K (m/s) err_K (m/s) planet mass (M_earth) error (M_earth) R (R_earth) #vs #phase1 #phase2 #phase3 #phase4 #phase5
        outfile.write(ostr)
        if writefit:
            outfp.write(ostr)


    if writefit:
        outfp.close()
    return



snames,gdfn,mfn,veldir,outdir,outfile,binned = parse_options()

gd = ascii.read(gdfn)
TESSAPFdata=ascii.read(mfn,format='csv')
for sname in snames:

    bvfn = "%sbinned.vels" % (sname)
    bsysfn = "%sbinned.sys" % (sname)    

    if binned == False:
        sfn = sname + ".sys"
        vfn = sname + ".vels"
        invels = fit_TESS_APF.readin_velsfile(os.path.join(veldir,vfn))
        ddates,dphases, dvels, derrs, di2sums = bin_phase_dates(invels["jd"],invels["phases"],invels['velocity'],invels["int unc"],invels["I2 counts"])
        if os.path.exists(os.path.join(veldir,bvfn)):
            os.remove(os.path.join(veldir,bvfn))
        ascii.write([ddates,dvels,derrs,di2sums,dphases], os.path.join(veldir,bvfn),format="no_header")

    
    indata = fit_TESS_APF.readin_velsfile(os.path.join(veldir,bvfn))
    planets, = np.where((TESSAPFdata['star_names'] == sname) & (TESSAPFdata['detected'] == "TRUE"))
    current_nvels =0
    for i in range(0,12):

        deltas = (indata['jd'] - indata['jd'].min() )< 90*(i+1)
        if len(indata['jd'][deltas]) <= current_nvels + 5:
            continue

        current_nvels = len(indata['jd'][deltas])
        
        tbvfn = "%sbinned_t.vels" % (sname)
        tbsysfn = "%sbinned_t.sys" % (sname)
        tname = "%sbinned_t" % (sname)

        if os.path.exists(os.path.join(veldir,tbvfn)):
            os.remove(os.path.join(veldir,tbvfn))
        if os.path.exists(os.path.join(veldir,tbsysfn)):
            os.remove(os.path.join(veldir,tbsysfn))
        
        ascii.write([indata['jd'][deltas],indata['velocity'][deltas],indata['int unc'][deltas],indata['I2 counts'][deltas],indata['phases'][deltas]], os.path.join(veldir,tbvfn),format="no_header")
        Generate_Velocities.write_sys(TESSAPFdata,sname,velname=tname,outdir=veldir)

        k=SystPy.Kernel()
        k.setEpoch(JD0)
        ddir =veldir
        k.addDataFile(tbsysfn, directory=ddir)
        fit_TESS_APF.add_planets(k,TESSAPFdata,planets,gd['vel_offset'][gd['starname'] == sname])
        
        print ("%s: RMS of fit %f" % (sname,k.getRms()))

        fit_planets(k,tname,TESSAPFdata['phase'][planets],TESSAPFdata['vmag'][planets],TESSAPFdata['rplanet'][planets],TESSAPFdata['true_mass'][planets],TESSAPFdata['Index'][planets],gd['vel_offset'][gd['starname'] == sname],outfile,writefit=True,veldir=veldir,outdir=outdir)


        
