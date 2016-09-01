import optparse
import os
import glob
import re

outpath = os.path.join(os.curdir,"..","SimFiles")

parser = optparse.OptionParser()
parser.add_option("-d","--date",dest="date",default="")
parser.add_option("-c","--cleanup",dest="cleanup",default=False,action="store_true")
(options, args) = parser.parse_args()    
if options.date != "":
    simfiles = glob.glob(options.date + ".simout")
else:
    simfiles = glob.glob("????-??-??.simout")

if options.cleanup:
    oldsims = glob.glob(outpath + "/*.sim")
    for osim in oldsims:
        os.unlink(osim)
    
ofp = False
for simfile in simfiles:
    simfp = open(simfile)
    lastline = ""
    for line in simfp:
        if re.match("\A\#",line):
            continue
        simdata = line.split()
        starname = simdata[0]
        if lastline != starname:
            if ofp:
                ofp.close()
            ofn = "%s.sim" % (starname)
            ofn = os.path.join(outpath,ofn)
            if os.path.exists(ofn):
                ofp = open(ofn,"a+")
            else:
                ofp = open(ofn,"w")
        outstr = "%s %s %s %s\n" % (simdata[3],simdata[5],simdata[6],simdata[7])
        ofp.write(outstr)
        lastline = starname
