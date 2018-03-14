import numpy as np
from astropy.io import ascii
from astropy.table import Table, Column

fn = "newgoogledex_R1.csv"
ofn = "newgoogledex_sinnoise.csv"

gd = ascii.read(fn)
gd['shortphase']=np.random.uniform(size=len(gd['vmag']))
start = 2.5 / (24.*60.)
end = 12.5 / (24.*60.)
gd['shortperiod']=np.random.uniform(start,end,size=len(gd['vmag']))


gd['longphase']=np.random.uniform(size=len(gd['vmag']))
gd['longperiod']=gd['rot_period']

gd['pri_offset']=np.zeros_like(gd['vmag'])
gd['vel_offset']=np.zeros_like(gd['vmag'])
gd['lastobs']=np.zeros_like(gd['vmag'])
gd['cadence']=np.zeros_like(gd['vmag'])
gd['cadence'] += 0.7

ascii.write(gd, output=ofn, format="csv")
