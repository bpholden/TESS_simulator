import optparse
import getpriority
from glob import glob
import re
import numpy as np
from astropy.io import ascii

def bin_phase_dates(dates,phases):

    dailydates = []
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
            print "crossed a phase boundary"
        curdates = dates[curday]
        dailydates.append(np.average(curdates))
        dailyphases.append(np.average(curphases))
        done[curday] = True
    return dailydates,dailyphases


parser = optparse.OptionParser()
parser.add_option("-p","--period",dest="period",default=0)
#parser.add_option("-s","--star",dest="starname",default="")
parser.add_option("-i","--infile",dest="infile",default="newgoogledex.csv")

(options, args) = parser.parse_args()

whole_filelist=glob('../SimFiles/*.sim')

googledex = ascii.read(options.infile)

phase_edges = [ getpriority.EDGE1, getpriority.EDGE2, getpriority.EDGE3, getpriority.EDGE4, getpriority.EDGE5]

for fn in whole_filelist:
    m = re.search("(TESSAPF\d+)\.sim",fn)
    if m:
        objn = m.group(1)
        match = googledex['starname'] == objn
        simvals = ascii.read(fn,names=["JD","i2","unc","dev"])
        phases = getpriority.compute_currentphase(simvals['JD'],googledex['foldperiod'][match],googledex['initialphase'][match])
        dailydates,dailyphases = bin_phase_dates(simvals['JD'],phases)
        print objn, dailydates,dailyphases
        binnum = np.digitize(dailyphases, phase_edges)
        for n in range(0,5):
            nobs = len(binnum[binnum == n])
            strn = "phase%dbin" % (n)
            googledex[strn][match] = nobs

ascii.write(googledex,options.infile,delimiter=",")
