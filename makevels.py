from astropy.io import ascii
import glob
from Generate_Velocities import readin_data
from Generate_Velocities import calc_rvs
from Generate_Velocities import write_vels
from Generate_Velocities import write_sys

#Read in TESSAPF master data file
file='../Datafiles/TESSAPF_AWMasses_prec_phase.csv'
TESSAPFdata=ascii.read(file,format='csv')

#Make a list of all of the APF simulation files we have
whole_filelist=glob.glob('../SimFiles/*.sim')

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
    write_vels(starname,simdata,velocities,phases)

    #write .sys file for each star
    write_sys(TESSAPFdata,starname)
