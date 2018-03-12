from __future__ import print_function

import numpy as np
import scipy as sp
from astropy.io import ascii
from astropy.table import Table, Column
from consts import JD0, ECC, OMEGA
import sim_nights

#Set number of observations required for a phase bin to be considered saturated
REQUIREDOBS=10

#Set edges for priority bins
EDGE0=0.0
EDGE1=.125
EDGE2=.375
EDGE3=.625
EDGE4=.875
EDGE5=1.


def compute_currentphase(currentJD,obsperiod,initialphase):
    bad = initialphase < 0
    currentphase=np.mod((np.mod(currentJD-JD0,obsperiod)/obsperiod)+initialphase,1)
    if bad.any():
        currentphase[bad] = -1

    return currentphase


def getpriority_inquad(starlist,data,currentJD,standard=False):
    """Takes in a list of strings (starlist), an astropy Table with the correct columns, and a specific JD (currentJD), then computes the priorities and RV phases of all the stars in starlist, returning them at the end of the function.
    starlist - a list of strings or an astropy Table Column of strings
    data - an astropy Table with all of the needed columns, frequently known as the googledex (ahem)
    currentJD - the Julian Date as a float that the phases and priorities will be computed for. 
    """
#
####################################################################################


    nstars=len(data)

    #Define new priority array that will contain each star's priority at t=JD
    priority=np.zeros(nstars)

    #Figure out where in the 0:1 RV phase curve each star is based on the time since JD0 and the period of the planet that is being used to set the priority
    currentphase = compute_currentphase(currentJD,data['foldperiod'],data['initialphase'])
    if standard:
        priority[data['initialphase']== -1] = 10.
        return priority,currentphase

    #Check the observation density of the stars' current phase bin
#    phasebin0sat=(data['phase0bin']/REQUIREDOBS)>=1
    phasebin1sat=(data['phase1bin']/REQUIREDOBS)>=1
    phasebin2sat=(data['phase2bin']/REQUIREDOBS)>=1
    phasebin3sat=(data['phase3bin']/REQUIREDOBS)>=1
    phasebin4sat=(data['phase4bin']/REQUIREDOBS)>=1

    quadsat = phasebin1sat + phasebin3sat

    binedges = [EDGE1,EDGE2,EDGE3,EDGE4,EDGE5]
    phasebins = np.digitize(currentphase,binedges,right=True)
            
    # this sets the priority for the output
    # this uses two sets of numbers to set the priority
    # first number is the current phase bin, the variable phasebins is based on only the current date
    # the second number is the number of past observations in the current phase bin, we only care about the number
    # divided by the number of desired observations 
#    priority[((phasebins == 1) & (phasebin1sat == 0)) | ((phasebins == 3) & (phasebin3sat == 0))] = 10
#    priority[((phasebins == 1) & (phasebin1sat == 1)) | ((phasebins == 3) & (phasebin3sat == 1))] = 3
#    priority[((phasebins == 0) & (phasebin0sat == 1)) | ((phasebins == 2) & (phasebin2sat == 1)) | ((phasebins == 4) & (phasebin4sat == 1))] = 1
#    priority[((phasebins == 0) & (phasebin0sat == 0)) | ((phasebins == 2) & (phasebin2sat == 0)) | ((phasebins == 4) & (phasebin4sat == 0))] = 5
#    priority[((phasebins == 0) & (phasebin0sat == 0) & (quadsat ==2)) | ((phasebins == 2) & (phasebin2sat == 0) & (quadsat ==2)) | ((phasebins == 4) & (phasebin4sat == 0) & (quadsat ==2))] = 8

    # remap priorities to a much simpler version
    priority[(phasebins == 1) | (phasebins == 3)] = 10
    priority[(phasebins == 0) | (phasebins == 2) | (phasebins == 4)] = 5
        
    #Case of standard star
    priority[data['initialphase']== -1] = 7
        
#    reallyinphase = ((currentphase > FIRSTSTART) & (currentphase < FIRSTEND)) | ((currentphase > SECONDSTART) & (currentphase < SECONDEND)) 
    
    #Add 1 to priorities if planet is <4 R_earth 
    priority[data['minradius'] < 4] += 1
#    priority[data['foldperiod'] > 30 ] += 1
#    priority[reallyinphase] += 1
#    priority += data['pri_offset']

    return priority,currentphase

def getpriority_outquad(starlist,data,currentJD,standard=False):
    """Takes in a list of strings (starlist), an astropy Table with the correct columns, and a specific JD (currentJD), then computes the priorities and RV phases of all the stars in starlist, returning them at the end of the function.
    starlist - a list of strings or an astropy Table Column of strings
    data - an astropy Table with all of the needed columns, frequently known as the googledex (ahem)
    currentJD - the Julian Date as a float that the phases and priorities will be computed for. 
    """
#
####################################################################################


    nstars=len(data)

    #Define new priority array that will contain each star's priority at t=JD
    priority=np.zeros(nstars)

    #Figure out where in the 0:1 RV phase curve each star is based on the time since JD0 and the period of the planet that is being used to set the priority
    currentphase = compute_currentphase(currentJD,data['foldperiod'],data['initialphase'])
    if standard:
        priority[data['initialphase']== -1] = 10.
        return priority,currentphase

    #Check the observation density of the stars' current phase bin
#    phasebin0sat=(data['phase0bin']/REQUIREDOBS)>=1
    phasebin1sat=(data['phase1bin']/REQUIREDOBS)>=1
    phasebin2sat=(data['phase2bin']/REQUIREDOBS)>=1
    phasebin3sat=(data['phase3bin']/REQUIREDOBS)>=1
    phasebin4sat=(data['phase4bin']/REQUIREDOBS)>=1

    quadsat = phasebin2sat + phasebin4sat

    binedges = [EDGE1,EDGE2,EDGE3,EDGE4,EDGE5]
    phasebins = np.digitize(currentphase,binedges,right=True)
            
    # this sets the priority for the output
    # this uses two sets of numbers to set the priority
    # first number is the current phase bin, the variable phasebins is based on only the current date
    # the second number is the number of past observations in the current phase bin, we only care about the number
    # divided by the number of desired observations 
    # priority[((phasebins == 1) & (phasebin1sat == 0)) | ((phasebins == 3) & (phasebin3sat == 0))] = 5
    # priority[((phasebins == 1) & (phasebin1sat == 1)) | ((phasebins == 3) & (phasebin3sat == 1))] = 1
    # priority[((phasebins == 0) & (phasebin0sat == 1)) | ((phasebins == 2) & (phasebin2sat == 1)) | ((phasebins == 4) & (phasebin4sat == 1))] = 3
    # priority[((phasebins == 0) & (phasebin0sat == 0)) | ((phasebins == 2) & (phasebin2sat == 0)) | ((phasebins == 4) & (phasebin4sat == 0))] = 10
    # priority[((phasebins == 1) & (phasebin1sat == 0) & (quadsat ==3)) | ((phasebins == 3) & (phasebin3sat == 0) & (quadsat ==3)) ] = 8

    priority[((phasebins == 1) ) | ((phasebins == 3) )] = 5
    priority[((phasebins == 0) ) | ((phasebins == 2) ) | ((phasebins == 4))] = 10
    
    #Case of standard star
    priority[data['initialphase']== -1] = 7
        
    
    #Add 1 to priorities if planet is <4 R_earth or if period > 30d
    priority[data['minradius'] < 4] += 1
#    priority[data['foldperiod'] > 30 ] += 1
#    priority += data['pri_offset']

    return priority,currentphase

def getpriority_uniform_old(starlist,data,currentJD,standard=False):
    """Takes in a list of strings (starlist), an astropy Table with the correct columns, and a specific JD (currentJD), then computes the priorities and RV phases of all the stars in starlist, returning them at the end of the function.
    starlist - a list of strings or an astropy Table Column of strings
    data - an astropy Table with all of the needed columns, frequently known as the googledex (ahem)
    currentJD - the Julian Date as a float that the phases and priorities will be computed for. 
    """
#
####################################################################################


    nstars=len(data)

    #Define new priority array that will contain each star's priority at t=JD
    priority=np.zeros(nstars)

    #Figure out where in the 0:1 RV phase curve each star is based on the time since JD0 and the period of the planet that is being used to set the priority
    currentphase = compute_currentphase(currentJD,data['foldperiod'],data['initialphase'])
    if standard:
        priority[data['initialphase']== -1] = 10.
        return priority,currentphase

    #Check the observation density of the stars' current phase bin
#    phasebin0sat=(data['phase0bin']/REQUIREDOBS)>=1
    phasebin1sat=(data['phase1bin']/REQUIREDOBS)>=1
    phasebin2sat=(data['phase2bin']/REQUIREDOBS)>=1
    phasebin3sat=(data['phase3bin']/REQUIREDOBS)>=1
    phasebin4sat=(data['phase4bin']/REQUIREDOBS)>=1

    binedges = [EDGE1,EDGE2,EDGE3,EDGE4,EDGE5]
    phasebins = np.digitize(currentphase,binedges,right=True)

    for i in range(0,5):
        namebin = "phase%dbin" %(i)
        priority[phasebins == i] = REQUIREDOBS - data[namebin][phasebins == i]
   
    #Case of standard star
    priority[data['initialphase']== -1] = 7
        
    
    #Add 1 to priorities if planet is <4 R_earth or if period > 30d
    priority[data['minradius'] < 4] += 1
#    priority[data['foldperiod'] > 30 ] += 1
#    priority += data['pri_offset']

    return priority,currentphase

def getpriority_uniform(starlist,data,currentJD,star_dates,standard=False):
    """Takes in a list of strings (starlist), an astropy Table with the correct columns, and a specific JD (currentJD), then computes the priorities and RV phases of all the stars in starlist, returning them at the end of the function.
    starlist - a list of strings or an astropy Table Column of strings
    data - an astropy Table with all of the needed columns, frequently known as the googledex (ahem)
    currentJD - the Julian Date as a float that the phases and priorities will be computed for. 
    """
#
####################################################################################


    nstars=len(starlist)

    #Define new priority array that will contain each star's priority at t=JD
    priority=np.zeros(nstars)

    #Figure out where in the 0:1 RV phase curve each star is based on the time since JD0 and the period of the planet that is being used to set the priority
    currentphase = compute_currentphase(currentJD,data['foldperiod'],data['initialphase'])
    #Case of standard star
    if standard:
        priority[data['initialphase']== -1] = 10.
        return priority,currentphase
    

    for i in range(0,nstars):
        if starlist[i] not in star_dates.keys():
            priority[i] = 10.
        else:
            if data['initialphase'][data['starname'] == starlist[i]] > 0:
                pastphases = compute_currentphase(np.asarray(star_dates[starlist[i]]),data['foldperiod'][data['starname'] == starlist[i]],data['initialphase'][data['starname'] == starlist[i]])
                pastphases = np.sort(pastphases)
                pastphases = np.append(pastphases, [pastphases[0]+1,pastphases[-1] - 1]) # wrap the end points, because we want to wrap the phase calculation
                # ie, if the current phase is 0.95 and then the closest could be 0.05
                pastphases -= currentphase[i]
                priority[i] = np.min(np.abs(pastphases)) * 10.0
            else:
                priority[i] = 1
    
    #Add 1 to priorities if planet is <4 R_earth or if period > 30d
    priority[data['minradius'] < 4] += 1
#    priority[data['foldperiod'] > 30 ] += 1
#    priority += data['pri_offset']

    return priority,currentphase

def getpriority_random(starlist,data,currentJD,standard=False):
    """Takes in a list of strings (starlist), an astropy Table with the correct columns, and a specific JD (currentJD), then computes the priorities and RV phases of all the stars in starlist, returning them at the end of the function.
    starlist - a list of strings or an astropy Table Column of strings
    data - an astropy Table with all of the needed columns, frequently known as the googledex (ahem)
    currentJD - the Julian Date as a float that the phases and priorities will be computed for. 
    """
#
####################################################################################


    nstars=len(data)

    #Define new priority array that will contain each star's priority at t=JD
    priority=np.float64(np.random.permutation(nstars))

    #Figure out where in the 0:1 RV phase curve each star is based on the time since JD0 and the period of the planet that is being used to set the priority
    currentphase = compute_currentphase(currentJD,data['foldperiod'],data['initialphase'])
    if standard:
        priority *= 0.0
        priority[data['initialphase']== -1] = 10.
        return priority,currentphase
    else:
        priority[data['initialphase']== -1] = 0.
 
    return priority,currentphase


def getpriority_hour(starlist,data,currentJD,lst,standard=False):
    """Takes in a list of strings (starlist), an astropy Table with the correct columns, and a specific JD (currentJD), then computes the priorities and RV phases of all the stars in starlist, returning them at the end of the function.
    starlist - a list of strings or an astropy Table Column of strings
    data - an astropy Table with all of the needed columns, frequently known as the googledex (ahem)
    currentJD - the Julian Date as a float that the phases and priorities will be computed for. 
    """
#
####################################################################################


    nstars=len(data)

    hourangle = np.abs(lst - (data['ra']*180./np.pi)/15.)
    hourangle = np.abs(hourangle)
    args = np.argsort(hourangle)[::-1]
    #Define new priority array that will contain each star's priority at t=JD
    ps=np.arange(0,len(data))
    priority = np.zeros(len(data),dtype=int)
    priority[args] += ps
    #Figure out where in the 0:1 RV phase curve each star is based on the time since JD0 and the period of the planet that is being used to set the priority
    currentphase = compute_currentphase(currentJD,data['foldperiod'],data['initialphase'])
    if standard:
        priority *= 0
        priority[data['initialphase']== -1] = 10.
        return priority,currentphase
    else:
        priority[data['initialphase']== -1] = 0
 
    return priority,currentphase


def getpriority(starlist,data,currentJD,star_dates,lst,standard=False,method="inquad"):
    if method == "uniform":
        priority,currentphase = getpriority_uniform(starlist,data,currentJD,star_dates,standard=standard)
        return priority,currentphase
    if method == "outquad":
        priority,currentphase = getpriority_outquad(starlist,data,currentJD,standard=standard)
        return priority,currentphase
    if method == "random":
        priority,currentphase = getpriority_random(starlist,data,currentJD,standard=standard)
        return priority,currentphase
    if method == "hour":
        priority,currentphase = getpriority_hour(starlist,data,currentJD,lst,standard=standard)
        return priority,currentphase

    priority,currentphase = getpriority_inquad(starlist,data,currentJD,standard=standard)
    return priority,currentphase


    

if __name__ == "__main__":
    star_strs = dict()
    star_dates = dict()
    mastername = "sim_master.simout"
    masterfp = open(mastername)
    for ln in masterfp:
        sim_nights.sim_results(ln,star_strs,star_dates)
    masterfp.close()
        
    file='../Datafiles/newgoogledex.csv'
    data=ascii.read(file,format='csv')
    currentJD = 2457470.285449536
    print ("in quad")
    print (getpriority(data['starname'],data,currentJD,star_dates,method="inquad"))
    print ("out of quad")
    print (getpriority(data['starname'],data,currentJD,star_dates,method="outquad"))
    print ("random")
    print (getpriority(data['starname'],data,currentJD,star_dates,method="random"))
    print ("uniform")
    print (getpriority(data['starname'],data,currentJD,star_dates,method="uniform"))

    print ("in quad")
    pri,cur= getpriority(data['starname'],data,currentJD,star_dates,standard=True,method="inquad")
    print (pri[pri > 0])
    print ("out of quad")
    pri,cur= getpriority(data['starname'],data,currentJD,star_dates,standard=True,method="outquad")
    print (pri[pri > 0])
    print ("random")
    pri,cur= getpriority(data['starname'],data,currentJD,star_dates,standard=True,method="random")
    print (pri[pri > 0]    )
    print ("uniform")
    pri,cur= getpriority(data['starname'],data,currentJD,star_dates,standard=True,method="uniform")
    print (pri[pri > 0]    )

