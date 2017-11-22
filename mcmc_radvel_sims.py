import os
import glob
import subprocess
import re

dns = glob.glob("../VelsFiles/36months*")

for dn in dns:
    veldir = dn
    pfdir = re.sub("VelsFiles","PlanetFitting",dn)
    exstr = "python fit_APF_TESS_radvel.py -v %s -o %s -a" % (veldir,pfdir)
    out = subprocess.check_output(exstr,shell=True)
    print exstr
    print out
