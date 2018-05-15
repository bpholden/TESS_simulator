from astropy.io import ascii
import glob
from Generate_Velocities import readin_data
from Generate_Velocities import calc_rvs
from Generate_Velocities import write_vels
from Generate_Velocities import write_sys
import optparse
import os

parser = optparse.OptionParser()
parser.add_option("-f","--file",dest="datafile",default="../Datafiles/TESSAPF_AWMasses_prec_phase_R1_vsinicut.csv")
parser.add_option("-i","--indir",dest="indir",default="../SimFiles/")        
parser.add_option("-o","--outdir",dest="outdir",default=".")        
(options, args) = parser.parse_args()    


if not os.path.isdir(options.outdir):
    os.mkdir(options.outdir)
 
#Read in TESSAPF master data file

TESSAPFdata=ascii.read(options.datafile,format='csv')

#Make a list of all of the APF simulation files we have
if options.indir[-1] != "/":
    options.indir =  options.indir + "/"
whole_filelist=glob.glob(options.indir + '*.sim')

#Determine which stars are TESSAPF Stars and need .vels and .sys files
sub='TESSAPF'
filelist=[s for s in whole_filelist if sub in s]

#For each star with a simulation file:
for i in range(0,len(filelist)):

    #read in simulation results
    simdata=readin_data(filelist[i])

    #pull out the star name
    temp=filelist[i].split('/')
    ll=len(temp)
    startemp=temp[ll-1].split('.')
    starname=startemp[0]

    #calculate the RV value for each observation in the .sim file
    velocities,phases=calc_rvs(TESSAPFdata,simdata,starname)
    
    #write .vels file for each star
    write_vels(starname,simdata,velocities,phases,outdir=options.outdir)

    #write .sys file for each star
    write_sys(TESSAPFdata,starname,outdir=options.outdir)
