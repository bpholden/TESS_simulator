from __future__ import print_function
from astropy.io import ascii
import os
import os.path
from glob import glob
import re
import optparse
import sys

import Generate_Velocities

def parse_fit(infile):
    masses = []
    errs = []
    ks = []
    errks = []
    rplanet = []
    nobs = []
    periods = []
    try:
        infp = open(infile)
    except:
        print ("cannot open %s" % (infile))
        return
    for ln in infp:
        d = ln.split()
        if len(d) < 7:
            continue
        else:
            periods.append(d[0])
            ks.append(d[1])
            errks.append(d[2])
            masses.append(d[3])
            errs.append(d[4])
            rplanet.append(d[5])
            nobs.append(d[6])
            
    return masses, errs, rplanet, nobs, periods, ks, errks

def findfiles(indir):
    bfits = glob(indir + "/" + "TESSAPF*.radvelfit")
    return bfits


def compile_data(indir,outfp,TESSAPFdata):
    fnames = findfiles(indir)
    sn = []
    ntot = 0
    nsn = 0
    nfifty = 0
    for fn in fnames:
        m = re.search("(TESSAPF\d+)",fn)
        sname = m.group(1)
        truemass = TESSAPFdata['true_mass'][(TESSAPFdata['star_names'] == sname) & (TESSAPFdata['detected']=='TRUE')]
        masses, errs, rplanets, nobs, periods, ks, errks = parse_fit(fn)
        np = len(masses)
        if int(nobs[0]) > 4:
            ntot += 1
        for i in range(0,np):
            ostr = "%s.%d " % (sname,(i+1))

            kmass = Generate_Velocities.calc_mass(ks[i],TESSAPFdata,sname)
            kerr = kmass * errks[i] / ks[i] 
            ostr += "%s %s %f %s %s %s %s %s\n" % (kmass,kerr,truemass[i],rplanets[i],periods[i],nobs[i],ks[i],errks[i])
            outfp.write(ostr)
            cmass = float(masses[i])
            csn = cmass / float(errs[i])
            sn.append(csn)
            if csn > 7.0:
                nsn += 1

            if cmass < truemass[i] + 0.5*truemass[i] and cmass > truemass[i] - 0.5*truemass[i]:
                nfifty += 1
    return sn, ntot, nsn,nfifty


parser = optparse.OptionParser()
parser.add_option("-i","--indir",dest="indir",default=".")
parser.add_option("-o","--outfile",dest="outfile",default="mass_radius")
parser.add_option("-t","--tessapf",dest="mfn",default='../Datafiles/TESSAPF_AWMasses_prec_phase_R1_vsinicut.csv')

(options, args) = parser.parse_args()
try:
    outfp = open(options.outfile,"w")
    hdr = "# starid mass(earths) err(earths) true_mass(earths) radius(earths) period(days) n_obs K(m/s) err(m/s) \n"
    outfp.write(hdr)
    outfp.write("#" + options.indir+"\n")
except Exception, e:
    print ("cannot open %s: %s" % (options.outfile,e))
    sys.exit(1)

try:
    TESSAPFdata=ascii.read(options.mfn,format='csv')
except Exception, e:
    print ("cannot open %s: %s" % (options.mfn, e))
    sys.exit(1)
    
sn, ntot, nsn, nfifty = compile_data(options.indir,outfp,TESSAPFdata)

print ("%d observed (n > 4 obs)" %(ntot))
print ("%d detected (sn > 7)" %(nsn))
print ("%d within 50 percent of true mass" %(nfifty))

