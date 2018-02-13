import numpy as np
from astropy.io import ascii
from astropy.table import Table, Column

fn = "newgoogledex.csv"
ofn = "newgoogledex_sinnoise.csv"

gd = ascii.read(fn)
gd['longphase']=np.random.uniform(size=len(gd['vmag']))
gd['shortphase']=np.random.uniform(size=len(gd['vmag']))
start = 2.5 / (24.*60.)
end = 12.5 / (24.*60.)
gd['shortperiod']=np.random.uniform(start,end,size=len(gd['vmag']))
start = 15.
end = 45.
gd['longperiod']=np.random.uniform(start,end,size=len(gd['vmag']))
ascii.write(gd, output=ofn, format="csv")
