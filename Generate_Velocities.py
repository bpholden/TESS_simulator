from __future__ import print_function
#from __future__ import division
import numpy as np
from astropy.io import ascii
from astropy.table import Table, Column
import scipy.constants as sc
import re
import os
from consts import JD0, ECC, OMEGA

def readin_data(infile):

    colnames = ["jd","i2cts","intunc","totunc"]

    try:
        simdata = ascii.read(infile,names=colnames)
    except:
        print("cannot open %s for input" % (infile))
        sys.exit(-1)

     
    return simdata

def calc_K(TESSAPFdata,starname,dotrue=True):
    
    if dotrue:
        loc=np.where(TESSAPFdata['star_names'] == starname)
        truemasses=TESSAPFdata['true_mass'][loc]
        #Pull out the planetary system info from the TESSAPF table
    else:
        loc=np.where((TESSAPFdata['star_names'] == starname) & (TESSAPFdata['detected'] == "TRUE"))
        truemasses=TESSAPFdata['est_mass'][loc]
        
    periods=TESSAPFdata['period'][loc]
    initphases=TESSAPFdata['phase'][loc]
    cosi=TESSAPFdata['cosi'][loc][0]    #same for every planet
    sini=np.sqrt(1-(cosi**2.))    #same for every planet
    mstar=TESSAPFdata['mstar'][loc][0]  #same for every planet (duh)

    #Convert to mks units
    period_sec=periods*86400.                   #days -> seconds
    truemasses=truemasses*5.97237*(10.**24.)  #Earth mass -> kg
    mstar= mstar*1988500*(10.**24.)           #Solar mass -> kg

    #Calulate the RV amplitude (K) of each planet
    K=((2*sc.pi*sc.G/period_sec)**(1./3.)) * (truemasses/((mstar+truemasses)**(2./3.)))

    return K

def calc_predrvs(TESSAPFdata,simdates,starname,dotrue=True):

    if dotrue:
        loc=np.where(TESSAPFdata['star_names'] == starname)
        truemasses=TESSAPFdata['true_mass'][loc]
        #Pull out the planetary system info from the TESSAPF table
    else:
        loc=np.where((TESSAPFdata['star_names'] == starname) & (TESSAPFdata['detected'] == "TRUE"))
        truemasses=TESSAPFdata['est_mass'][loc]
        
    periods=TESSAPFdata['period'][loc]
    initphases=TESSAPFdata['phase'][loc]
    cosi=TESSAPFdata['cosi'][loc][0]    #same for every planet
    sini=np.sqrt(1-(cosi**2.))    #same for every planet
    mstar=TESSAPFdata['mstar'][loc][0]  #same for every planet (duh)

    #Convert to mks units
    period_sec=periods*86400.                   #days -> seconds
    truemasses=truemasses*5.97237*(10.**24.)  #Earth mass -> kg
    mstar= mstar*1988500*(10.**24.)           #Solar mass -> kg

    #Calulate the RV amplitude (K) of each planet
    K=((2*sc.pi*sc.G/period_sec)**(1./3.)) * (truemasses/((mstar+truemasses)**(2./3.)))

    #Calculate the true anomaly (f) for each planet, which is == to the eccentric anamoly (u) here as we have set ecc=0

    currentphases=np.empty((len(simdates),len(periods)))
    for j in range(0,len(periods)):        
        currentphases[:,j]=np.mod((np.mod(simdates-JD0, periods[j])/periods[j])+initphases[j],1)

    f=currentphases*2.*sc.pi

    #Calculate the instantaneous radial velocity for each planet
    #Need to make this work for multiple planets
    Vinsts=np.empty_like(currentphases)
    for k in range(0,len(periods)):
        Vinsts[:,k]=K[k]*(np.cos(f[:,k]+OMEGA) + ECC*np.cos(OMEGA))

    Vtots=np.empty(len(currentphases))
    for l in range (0,len(currentphases)):
        Vtots[l]=Vinsts[l,:].sum()
    Vtots -= Vtots[0]
    return Vtots,currentphases



def calc_rvs(TESSAPFdata,simdata,starname,dotrue=True,phases=[]):
    

    if dotrue:
        #Pull out the planetary system info from the TESSAPF table
        loc=np.where(TESSAPFdata['star_names'] == starname)
        truemasses=TESSAPFdata['true_mass'][loc]
    else:
        #Pull out the planetary system info from the TESSAPF table
        loc=np.where((TESSAPFdata['star_names'] == starname) & (TESSAPFdata['detected'] == "TRUE"))
        truemasses=TESSAPFdata['est_mass'][loc]
        
    periods=TESSAPFdata['period'][loc]
    initphases=TESSAPFdata['phase'][loc]
    cosi=TESSAPFdata['cosi'][loc][0]    #same for every planet
    sini=np.sqrt(1-(cosi**2.))    #same for every planet
    mstar=TESSAPFdata['mstar'][loc][0]  #same for every planet (duh)

    #Convert to mks units
    period_sec=periods*86400.                   #days -> seconds
    truemasses=truemasses*5.97237*(10.**24.)  #Earth mass -> kg
    mstar= mstar*1988500*(10.**24.)           #Solar mass -> kg

    #Calulate the RV amplitude (K) of each planet
    K=((2*sc.pi*sc.G/period_sec)**(1./3.)) * (truemasses/((mstar+truemasses)**(2./3.)))

#    print('RV amplitude of planets:')
#    print(K)
#    print()

    #Calculate the true anomaly (f) for each planet, which is == to the eccentric anamoly (u) here as we have set ecc=0

    if len(phases) > 0:
        currentphases = phases
    else:
        currentphases=np.empty((len(simdata),len(periods)))
        for j in range(0,len(periods)):        
            currentphases[:,j]=np.mod((np.mod(simdata['jd']-JD0, periods[j])/periods[j])+initphases[j],1)

    
    f=currentphases*2.*sc.pi

#    print('true anomalies:')
#    print(f)
#    print()


    #Calculate the instantaneous radial velocity for each planet
    #Need to make this work for multiple planets
    Vinsts=np.empty((len(simdata),len(periods)))
    for k in range(0,len(periods)):
        Vinsts[:,k]=K[k]*(np.cos(f[:,k]+OMEGA) + ECC*np.cos(OMEGA))

#    print('instantaneous RV amplitudes:')
#    print(Vinsts)
#    print()

    #Calculate total radial velocity of the star (orbited by nplanets)
    Vtots=np.empty(len(simdata))
    for l in range (0,len(simdata)):
        Vtots[l]=Vinsts[l,:].sum()+simdata['totunc'][l]
    # RVs are relative
    Vtots -= Vtots[0]
#    print('Total instantaneous RVs of the star at specified JD:')
#    print(Vtots)

    return Vtots,currentphases

def write_vels(starname,simdata,velocities,phases,outdir="../VelsFiles/"):
    filename=starname+'.vels'
    filename = os.path.join(outdir,filename)
    f = open(filename,'w')
    for k in range(0,len(simdata)):
        ph=phases[k,:].astype(str)
        strph=' '.join(ph)
        if velocities[k] < 1000.0 and simdata['intunc'][k] < 1000.0 and simdata['intunc'][k] > 0.0:
            outstr = "%s %s %s %s %s \n" % (simdata['jd'][k],velocities[k],simdata['intunc'][k],simdata['i2cts'][k],strph)
        f.write(outstr)
    f.close()

    print ('vels file saved as: '+filename)
    
def write_sys(TESSAPFdata,starname,velname="",outdir='../VelsFiles/'):
    loc=np.where(TESSAPFdata['star_names'] == starname)
    starmass=TESSAPFdata['mstar'][loc][0]

    if velname != "":
        filename=velname+'.sys'
    else:
        filename=starname+'.sys'
    f = open(os.path.join(outdir,filename),'w')
    f.write('Data { \n')
    if velname != "":
        f.write('    RV[] "'+velname+'.vels" \n' )
    else:
        f.write('    RV[] "'+starname+'.vels" \n' )
            
    f.write('} \n')
    f.write('"'+starname+'" { \n')
    f.write(' Mass '+starmass.astype(str)+' \n')
    f.write('}')
    f.close()

    print ('sys file saved as: '+filename)
    
