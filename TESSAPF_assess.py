from __future__ import print_function
from __future__ import division
import numpy as np 
import scipy as sp
import matplotlib
from astropy.io import ascii
from astropy.table import Table, Column
import scipy.constants as sc
import glob
import csv
import re

PUBRATIO=7.0
DROPRATIO=1.5

def all_same(items):
    return all(x == items[0] for x in items)

def assess(input_gd_file,new_gd_file,fitfile_loc):
    #readin the googledex used for this run
    gd=ascii.read(input_gd_file,format='csv')
    
    #Read in the fits files of from an observing simulation
    fit_filelist=glob.glob(fitfile_loc+'/*binned.fit')
    print(fit_filelist)
    #Get list of all the TESSAPF star numbers that have .fit files
    starnums=[]
    for s in fit_filelist:
        mtch = re.search('F(\d+)b',s)
        if mtch:
            val = mtch.group(1)
            # if you need it as an integer
            intval = int(val)  
            starnums.append(intval)
    starlist=np.unique(starnums)
    
    #Run through each planet individually and assess how its velocities are doing
    pl_assess=[]
    needs_review=[]
    assess_names=[]

    for i in range(0,len(starlist)):
        sname='TESSAPF'+starlist[i].astype(str)
        file=glob.glob(fitfile_loc+'/'+sname+'binned.fit')
        gd_loc = [q for q,x in enumerate(gd['starname']) if x == sname]

    #read in info on fit to the RV residuals after known
    #planet signals have been removed
        f=open(file[0],'r')
        residfit=f.readline().split(' ')
        residfit_obsspan=float(residfit[0])
        residfit_slope=float(residfit[1])
        residfit_intercept=float(residfit[2])
        residfit_rms=float(residfit[3])
        f.close()

        #read in info on the individual planet fits
        fitdata=np.loadtxt(file[0],delimiter=' ',ndmin=2,skiprows=1)
        nplanets=len(fitdata)

        for k in range(0,nplanets):
            period=fitdata[k,0].astype(float)
            semiamp=fitdata[k,1].astype(float)
            semiamp_err=fitdata[k,2].astype(float)
            plmass_me=fitdata[k,3].astype(float)
            plmass_err=fitdata[k,4].astype(float)
            plradii=fitdata[k,5].astype(float)
            nbinvels=fitdata[k,6].astype(float)
            nbinq1=fitdata[k,8].astype(float)
            nbinq2=fitdata[k,10].astype(float)

            fitratio=semiamp/semiamp_err

            #First check if there's evidence of a linear trend (long per compation)
            #If so, mark the planet as needing review and move on
            if (residfit_rms > 2*max(1,gd['precision'][gd_loc[0]])):
                stat='monitor'  #used to be needs_review
                pl_assess.append(['TESSAPF'+starlist[i].astype(str),period,semiamp,nbinvels,nbinq1+nbinq2,stat,plmass_me,plmass_err,plradii]) 
                assess_names.append('TESSAPF'+starlist[i].astype(str))

#                needs_review.append(sname) 
                print(['TESSAPF'+starlist[i].astype(str),period,stat])
                continue
            
            if (nbinq1 < 10) | (nbinq2 < 10) | ((fitratio >= DROPRATIO) & (fitratio < PUBRATIO)): 
                stat='monitor'
            if (fitratio < DROPRATIO) & (nbinq1 >= 10) & (nbinq2 >= 10): #was previously all binned vels
                stat='drop'
            if fitratio >= PUBRATIO: 
                stat='publish'
            pl_assess.append(['TESSAPF'+starlist[i].astype(str),period,semiamp,nbinvels,nbinq1+nbinq2,stat,plmass_me,plmass_err,plradii])
            assess_names.append('TESSAPF'+starlist[i].astype(str))
            print(['TESSAPF'+starlist[i].astype(str),period,stat,plmass_me,plmass_err])


    pri_offset=np.zeros(len(gd))
    new_periods=np.zeros(len(gd))
    new_Ks=np.zeros(len(gd))

    publishable=[]
    monitor=[]
    drop=[]
    needs_review=[]
    changed_fold=[]
    changed_prec=[]
    changed_priority=[]
    DropNbin=20

    for i in range (0,len(starlist)):
    
        sname='TESSAPF'+starlist[i].astype(str)
        gd_loc = [qq for qq,x in enumerate(gd['starname']) if x == sname]
        loc = [qq for qq,x in enumerate(assess_names) if x == sname]

        periods=[]
        Ks=[]
        checknqvels=[]
        stat=[]
    
        for j in range(0,len(loc)):
            assessment=pl_assess[loc[j]]
            periods.append(assessment[1])
            Ks.append(assessment[2])
            checknqvels.append(assessment[4])
            stat.append(assessment[5])


        #If any of the planets have a "need review" tag
        #add the star to the Review array and move to next star
        if 'needs_review' in stat:
            pri_offset[gd_loc]=0
            needs_review.append(sname)
            continue

        #Case of all planets in system having the same status:
        if (np.size(stat) == 1) or (all_same(stat)):
            currentprec=gd['precision'][gd_loc]

            #Publish: set -5 priority offset for stars with one planet 
            #that is deemed publishable. Set new K precision if desired prec
            #in the googledex isn't the same as 1/2 the smallest fitted K.
            if 'publish' in stat:
                if min(Ks)/2. != currentprec:
                    new_Ks[gd_loc]=max(1,min(Ks))
                pri_offset[gd_loc]=-5
                publishable.append(sname)
            
            #Monitor: set 0 priority offset (so no change)for stars with one planet
            #that requires further monitoring. Set new K precision if desired prec
            #in the googledex isn't the same as 1/2 the smallest fitted K.
            if 'monitor' in stat:
                if min(Ks)/2. != currentprec:
                    new_Ks[gd_loc]=max(1,min(Ks))
                pri_offset[gd_loc]=0
                monitor.append(sname)
            
            #Drop: set -10 priority offset for stars with one planet that 
            #show K/err(K) < 1 after 10+ binned observations
            if 'drop' in stat:
                pri_offset[gd_loc]=-10
                drop.append(sname)
            

    #Case of all planets in system not having the same status:
        else:
        
            #Make note if the star has at least 1 publishable planet:
            if 'publish' in stat:
                publishable.append(sname)
        
            #Choose new period based on the planets that are neither publishable nor drop-worthy.
            #if there is such a planet, switch the fold period to its period and
            #add the star to the publish list and the changed_periods list
            monloc=[qq for qq,x in enumerate(stat) if x == 'monitor']
            monper = [periods[qq] for qq in monloc]

            pubmonloc=[i for i,x in enumerate(stat) if x != 'drop']
            pubmonper=[periods[qq] for qq in pubmonloc]
        
        
            #Set the desired precision to be 1 or the smallest of the monitored/publishable
            #planets suspected K values
            viableKs= [Ks[qq] for qq in pubmonloc] 
            new_Ks[gd_loc]=max(1,min(viableKs))
        
            if len(monloc) != 0:
            #If there are planets that need monitoring in the system, set the fold period
            #to the shortest period of these planets and set the priority offset to 0
                new_periods[gd_loc]=min(monper)
                pri_offset[gd_loc]=0
            
            else:
            #if there's not such a planet (everything is either publishable or should
            #be dropped) then set the priority offset to -5, and the period should be set 
            #to the shortest period of the publishable planets
                new_periods[gd_loc]=min(pubmonper)
                pri_offset[gd_loc]=-5
                drop.append(sname)
                publishable.append(sname)

            
    #update the observing precision and fold period values in the googledex
    #if this iteration already has the observing
    for i in range(0,len(gd)):
        if new_Ks[i] != 0:
            gd['precision'][i] = new_Ks[i] 
            changed_prec.append(gd['starname'][i])
        if new_periods[i] != 0:
            gd['foldperiod'][i]=new_periods[i]
            changed_fold.append(gd['starname'][i])
        if pri_offset[i] != 0:
            gd['pri_offset'][i]=pri_offset[i]
            changed_priority.append(gd['starname'][i])

    #write the new googledex to file
    ascii.write(gd,new_gd_file,format='csv')


    #Print the stars that fall into each of the status catagories
    print(len(publishable),' stars have at least one planet ready for publication:',publishable)
    print('')
    print(len(monitor),' stars need continued monitoring: ',monitor)
    print('')
    print(len(drop),' stars should be dropped from the program: ',drop)
    print('')
    print(len(needs_review),' stars need human review: ',needs_review)
    print('')
    print(len(changed_fold),' stars have a new fold period: ',changed_fold)
    print('')
    print(len(changed_prec),' stars have a new precision: ',changed_prec)
    print('')
    print(len(changed_priority),' stars have a new priority offset: ',changed_priority)


import optparse
import os
import sys
parser = optparse.OptionParser()
parser.add_option("-o","--outgoogledex",dest="outfile",default="newgoogledex_post.csv")
parser.add_option("-i","--ingoogledex",dest="infile",default="newgoogledex.csv")
parser.add_option("-p","--path",dest="path",default="../PlanetFitting/")
(options, args) = parser.parse_args()    

if not os.path.exists(options.infile):
    ostr= "%s does not exist" % (options.infile)
    print( ostr)
    sys.exit(-1)
    
if os.path.exists(options.outfile):
    print( "%s does exist, will not overwrite" % (options.outfile))
    sys.exit(-1)

if not os.path.isdir(options.path):
    print( "%s is not a directory" % (options.path))
    sys.exit(-1)

assess(options.infile,options.outfile,options.path)
