# Sept 1 2018

from datetime import datetime, timedelta
import random
import optparse
import subprocess

from sim_nights import gen_datelist

parser = optparse.OptionParser()
parser.add_option("-s","--start",dest="start",default="2019/8/30")
parser.add_option("-e","--end",dest="end",default="2019/9/1")
parser.add_option("-d","--double",dest="double",default=False,action="store_true")

(options, args) = parser.parse_args()    

dl = gen_datelist(options.start,options.end,double=options.double)
print dl
