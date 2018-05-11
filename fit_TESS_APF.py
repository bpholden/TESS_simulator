from __future__ import print_function
import re
import sys
import optparse
import os.path
from glob import glob
import re

from astropy.io import ascii
from astropy.table import Table
import numpy as np
import scipy as sp
import scipy.constants as sc
import scipy.optimize as op
import matplotlib.pyplot as plt

sys.path.append("../simulator/")
import getpriority
import Generate_Velocities
from consts import JD0, ECC, OMEGA
import SystPy

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
    parser.add_option("-i","--infile",dest="infile",default="../Datafiles/newgoogledex.csv")
    parser.add_option("-v","--veldir",dest="veldir",default="../VelsFiles/")
    parser.add_option("-o","--outdir",dest="outdir",default="../PlanetFitting/")        
    parser.add_option("-t","--tessapf",dest="mfn",default='../Datafiles/TESSAPF_AWMasses_prec_phase.csv')
    (options, args) = parser.parse_args()
    if len(args) < 1 and not options.all:
        print ("needs a star name")
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
                print ("%s does not exist" %(sname))
                sys.exit()
            snames.append(sname)

    if not os.path.isdir(options.outdir):
        try:
            os.mkdir(options.outdir)
        except Exception as e:
            print ("cannot make dir %s: %s" % (options.outdir,e))
            sys.exit()
            

    return snames,options.infile,options.mfn,veldir,options.outdir


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
        if curphases.max() - curphases.min() > 0.1:
            print ("crossed a phase boundary")
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



def add_planets(k,TESSAPFdata,planets,veloff):

    i = 1
    for idx in planets:
        ma = TESSAPFdata['phase'][idx] * 360.0 + OMEGA * 180. / sc.pi
        if ma < 0:
            ma = ma + 360.0
        k.addPlanet([SystPy.K_PER, TESSAPFdata['period'][idx], SystPy.K_MASS, TESSAPFdata['est_mass'][idx], SystPy.K_MA, ma, SystPy.K_ECC, 0.0, SystPy.K_DONE])

        for w in range(SystPy.K_ELEMENTS_SIZE):
            k.setElementFlag(i, w, SystPy.K_ACTIVE & SystPy.K_MINIMIZE)
        k.setElementFlag(i, SystPy.K_MASS, SystPy.K_ACTIVE | SystPy.K_MINIMIZE)
#        k.setElementFlag(i, SystPy.K_MA, SystPy.K_ACTIVE | SystPy.K_MINIMIZE)
        k.setPar(0,veloff) # sets the offset of the data set to be tabulated offset in m/s
        k.setParFlag(0, SystPy.K_MINIMIZE | SystPy.K_ACTIVE)       # sets it to be fit
    
        k.calculate()
        k.minimize()
        i+=1
    k.setTrendFlag(SystPy.K_ACTIVE | SystPy.K_MINIMIZE)
    k.calculate()
    k.minimize()
    return

def fit_dates_vels(dates,vels,errs):

    odates = dates - np.min(dates)
    try:
        pcoeff, cov = np.polyfit(odates,vels,1,w=1.0/errs**2,cov=True)
#    print (pcoeff)
#    txt = "%.4f %.4f %.4f %.4f %.4f" % (pcoeff[0],pcoeff[1],np.max(odates),np.sqrt(cov[0][0]),np.sqrt(cov[1][1]))
        txt = "%.4f %.4f %.4f" % (pcoeff[0],pcoeff[1],np.max(odates))
    except:
        txt = "Cannot perform polyfit"
    return txt


def checknumvels(invels,sname,outdir="../PlanetFitting"):
    kernelname = "%s" % (sname)
    kernelname = os.path.join(outdir,kernelname)
    if os.path.exists(kernelname):
        kernel = SystPy.loadKernel(kernelname)
        if kernel.nData() == len(invels['jd']):
            return False
    return True

def plot_planets(k,sname,initphases,vmag,rplanets,mplanets,planetids,veloff,writefit=False,veldir="../VelsFiles",outdir="../PlanetFitting"):
    
    fn = sname + ".vels"
    indata = readin_velsfile(os.path.join(veldir,fn))

    m = SystPy.matrix_to_array(SystPy.getResidualMatrix(k))
    txt = fit_dates_vels(m[:,0],m[:,1],m[:,2])
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
        daterange=np.linspace(np.min(indata['jd']),np.min(indata['jd'])+fperiod,num=200)
        plotting_phases=getpriority.compute_currentphase(daterange,fperiod,initphases[i])
        plotting_phases = plotting_phases[(plotting_phases > 0) & (plotting_phases < 1)]
        plotting_phases = np.sort(plotting_phases)
        
        plt.errorbar(phases,indata['velocity']-veloff,yerr=indata['int unc'],fmt='bo')
        
        plt.xlabel('Phases',fontsize=14)
        plt.ylabel('Velocity (m s$^{-1}$)',fontsize=14)
        plt.title('Phase folded velocities for ' + sname)

        #        K = k.getElement(i+1,SystPy.K_SEMIAMP)
        K = meds[i+1,SystPy.K_SEMIAMP]
        err_K = mads[i+1,SystPy.K_SEMIAMP]
        f=plotting_phases*2.*sc.pi
        pvels=(np.cos(f+OMEGA) + ECC*np.cos(OMEGA))
        plt.plot(plotting_phases,K*pvels,'k-')
        top = K*pvels+err_K
        bot = K*pvels-err_K
        plt.fill_between(plotting_phases,bot,top,color='blue',alpha=0.2)

        ylims = plt.ylim()
        plt.xlim(0,1)
        xs = np.asarray([getpriority.EDGE1,getpriority.EDGE2])
        plt.fill_between(xs,ylims[0],ylims[1],facecolor='grey',alpha=0.2)
        xs = np.asarray([getpriority.EDGE3,getpriority.EDGE4])
        plt.fill_between(xs,ylims[0],ylims[1],facecolor='grey',alpha=0.2)
        pstr= "Period %.2f days M=%.2f M$_{sun}$ V = %.2f mag" % (fperiod,k.getElement(0,1),vmag[i])
        plt.text(0.37,ylims[1]*0.9,pstr,va='bottom')
        pstr= "$M=%.2g \pm %.2g\ M_J$"% (k.getElement(i+1,SystPy.K_MASS),mads[i+1,SystPy.K_MASS])
        plt.text(0.37,ylims[1]*0.8,pstr,va='bottom')        
        pstr="$(K = %.2f \pm %.2f \ m\ s^{-1})$" %(K,err_K)
        plt.text(0.37,ylims[1]*0.7,pstr,va='bottom')
        pstr= "$\sigma$=%.2f m s$^{-1}$ $\mu$=%.2f deg" % (k.getRms(),k.getElement(i+1,SystPy.K_MA))
        plt.text(0.37,ylims[1]*0.6,pstr,va='bottom')

        binedges = [getpriority.EDGE1,getpriority.EDGE2,getpriority.EDGE3,getpriority.EDGE4,getpriority.EDGE5]
        phasebins = np.digitize(phases,binedges,right=True)
        ninbins = []
        for nbin in range(0,len(binedges)):
            cpb = phasebins[phasebins == nbin]
            ninbins.append(len(cpb))

        figname = "%s_planet%d_PhasedVelocities.pdf" % (sname,i+1)
        figname = os.path.join(outdir,figname)
        plt.savefig(figname, bbox_inches='tight')
        plt.close()
        if writefit:
            ostr = "%f %f %f %.4g %.4g %.4g %d" % (fperiod,K,err_K,k.getElement(i+1,SystPy.K_MASS)*317.83,mads[i+1,SystPy.K_MASS]*317.83,rplanets[i],len(phases))
            # period (days) K (m/s) err_K (m/s) planet mass (M_earth) error (M_earth) R (R_earth) #vs #phase1 #phase2 #phase3 #phase4 #phase5
            bstr = ""
            for i in range(0,len(ninbins)):
                bstr += " %d" % ninbins[i]
            bstr += "\n"
            outfp.write(ostr)
            outfp.write(bstr)

    kernelname = "%s" % (sname)
    kernelname = os.path.join(outdir,kernelname)
    k.save(kernelname)
    if writefit:
        outfp.close()
    return

def make_periodogram(k,sname,outdir="../PlanetFitting"):
    minPer=1.5 # ignore aliasing on short time scales
    p = SystPy.periodogram_ls(SystPy.getResidualMatrix(k), 100000, minPer, 10000, 0, SystPy.K_T_TIME, SystPy.K_T_SVAL, SystPy.K_T_ERR)
    amp = p[:,SystPy.K_PS_Z]
    peaks = SystPy.getPeakIndx(amp)
    amp = amp[peaks]
    per = p[:,SystPy.K_PS_TIME][peaks]
    fap = p[:,SystPy.K_PS_FAP][peaks]
    sort = amp.argsort()[::-1]
    #print (amp[sort])
    #print (per[sort])
    #print (fap[sort])
    t = Table([amp[sort],per[sort],fap[sort]], names=['Amp','Period','FAP'])
#    print (t[0:5])

    
    plt.plot(p[:,SystPy.K_PS_TIME], p[:,SystPy.K_PS_Z], c='black')
    plt.xscale('log')
    plt.xlabel("Period", fontsize=14)
    plt.ylabel("Power", fontsize=14)
    plt.xlim([0.5, 10000])
    ylims = plt.ylim()
    plt.title("Periodgram of residuals for %s" % (sname))
    ttxt = "Amp Period FAP" 
    plt.text(1000,0.95*ylims[1],ttxt)
    for i,s in enumerate([0.9,0.85,0.8]):
        ttxt = "%.2f %.2f %.2g" % (t['Amp'][i],t['Period'][i],t['FAP'][i])
        plt.text(1000,s*ylims[1],ttxt)
    pname = "%s_ResidualPeriodogram.pdf" % (sname)
    pname = os.path.join(outdir,pname)    
    plt.savefig(pname, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":

    snames,gdfn,mfn,veldir,outdir = parse_options()

    gd = ascii.read(gdfn)
    TESSAPFdata=ascii.read(mfn,format='csv')
    for sname in snames:

        sfn = sname + ".sys"
        vfn = sname + ".vels"
        newvels = True

        invels = readin_velsfile(os.path.join(veldir,vfn))
        ddates,dphases, dvels, derrs, di2sums = bin_phase_dates(invels["jd"],invels["phases"],invels['velocity'],invels["int unc"],invels["I2 counts"])
        bvfn = sname + "binned"
        ascii.write([ddates,dvels,derrs,di2sums,dphases], os.path.join(veldir,bvfn+".vels"),format="no_header",overwrite=True)
        Generate_Velocities.write_sys(TESSAPFdata,sname,velname=bvfn,outdir=veldir)

        newvels = checknumvels(invels,sname,outdir=outdir)
    
        if len(dvels) > 4 and newvels:
            k=SystPy.Kernel()
            k.setEpoch(JD0)
            ddir =veldir
            k.addDataFile(sfn, directory=ddir)

            planets, = np.where((TESSAPFdata['star_names'] == sname) & (TESSAPFdata['detected'] == "TRUE"))
            # planet masses are in earth masses, Systemic likes Jupiters
            #mearth = 5.9722e24
            #mjup = 1.898e27
            mratio = 317.83
            TESSAPFdata['est_mass'] /= mratio

            add_planets(k,TESSAPFdata,planets,gd['vel_offset'][gd['starname'] == sname])
    
            #print (k.getElements())
            print ("%s: RMS of fit %f" % (sname,k.getRms()))

            plot_planets(k,sname,TESSAPFdata['phase'][planets],TESSAPFdata['vmag'][planets],TESSAPFdata['rplanet'][planets],TESSAPFdata['true_mass'][planets],TESSAPFdata['Index'][planets],gd['vel_offset'][gd['starname'] == sname],veldir=veldir,outdir=outdir)
            # feature request, spit out periodogram
            make_periodogram(k,sname,outdir=outdir)
        # now do the binned velocities

            k=SystPy.Kernel()
            k.setEpoch(JD0)
            ddir =veldir
            bsfn = bvfn + ".sys"
            k.addDataFile(bsfn, directory=ddir)
            add_planets(k,TESSAPFdata,planets,gd['vel_offset'][gd['starname'] == sname])
    
            #print (k.getElements())
            print ("%s: RMS of fit %f" % (sname,k.getRms()))

            plot_planets(k,bvfn,TESSAPFdata['phase'][planets],TESSAPFdata['vmag'][planets],TESSAPFdata['rplanet'][planets],TESSAPFdata['true_mass'][planets],TESSAPFdata['Index'][planets],gd['vel_offset'][gd['starname'] == sname],writefit=True,veldir=veldir,outdir=outdir)
        
