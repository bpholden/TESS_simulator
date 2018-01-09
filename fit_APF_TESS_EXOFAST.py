from astropy.io import ascii
from astropy.table import Table
import numpy as np
import scipy as sp
import pandas as pd
import scipy.constants as sc
import scipy.optimize as op
import matplotlib.pyplot as plt
import re
import sys
import optparse
import copy
import os.path
from glob import glob
import re
sys.path.append("../simulator/")
import getpriority
import Generate_Velocities
from consts import JD0, ECC, OMEGA
import ConfigParser


def parse_options():

    parser = optparse.OptionParser()
    parser.add_option("-a","--all",dest="all",default=False,action="store_true")
    parser.add_option("-i","--infile",dest="infile",default="../Datafiles/newgoogledex.csv")
    parser.add_option("-v","--veldir",dest="veldir",default="../VelsFiles/")
    parser.add_option("-o","--outdir",dest="outdir",default="../PlanetFitting/")        
    parser.add_option("-t","--tessapf",dest="mfn",default='../Datafiles/TESSAPF_AWMasses_prec_phase.csv')
    parser.add_option("-e","--extra",dest="extra",default=False)
    (options, args) = parser.parse_args()
    if len(args) < 1 and not options.all:
        print "needs a star name"
        sys.exit()
    snames = []

    veldir = options.veldir + "/"
    
    if options.all:
        files = glob(veldir+"TESSAPF*.vels")
        snames = []
        for f in files:
            m = re.search("(TESSAPF\d+)\.vels",f)
            if m:
                snames.append(m.group(1))
    else:
        for a in args:
            if re.match("\A\d+\Z",a):
                sname = "TESSAPF%s" % (a)
            else:
                sname = a
            vname = os.path.join(veldir,sname + ".vels")
            if not os.path.exists(vname):
                print "%s does not exist" %(sname)
                sys.exit()
            snames.append(sname)

    if not os.path.isdir(options.outdir):
        try:
            os.mkdir(options.outdir)
        except Exception as e:
            print "cannot make dir %s: %s" % (options.outdir,e)
            sys.exit()
            

    return snames,options.infile,options.mfn,veldir,options.outdir,options.extra

def initialize_model(planets,TESSAPFdata,addextra=False):
    time_base = 2458362.5

    nplan = len(planets)
    if addextra:
        nplan = nplan + 1
    params = radvel.Parameters(nplan,basis='per tp e w k') # number of planets = nplan
    n = 1
    for idx in planets:
        c_period = TESSAPFdata['period'][idx]
        k = radvel.utils.semi_amplitude(TESSAPFdata['est_mass'][idx],c_period,TESSAPFdata['mstar'][idx],0.)
        tpstart = time_base + c_period * (TESSAPFdata['phase'][idx] + 0.25)
        params['per%d' % (n)] = radvel.Parameter(value=c_period)
        params['tp%d' % (n)] = radvel.Parameter(value=tpstart)
        params['e%d' % (n)] = radvel.Parameter(value=0.0)
        params['w%d' % (n)] = radvel.Parameter(value=0.)
        params['k%d' % (n)] = radvel.Parameter(value=k)
        n = n + 1

    if addextra:
        params['per%d' % (n)] = radvel.Parameter(value=30.)
        params['tp%d' % (n)] = radvel.Parameter(value=time_base + 15.)
        params['e%d' % (n)] = radvel.Parameter(value=0.0001)
        params['w%d' % (n)] = radvel.Parameter(value=np.pi/2.)
        params['k%d' % (n)] = radvel.Parameter(value=2.)
        
#     nparams = params.basis.from_cps(params,'per tc e w k') # transform to new basis (for less biased fits?)
#     n = 1
#     if addextra:
#         n = 2
#     for i in range(1,len(planets+n)):
# #         del nparams['e%d' % (i)]
# # #        del nparams['k%d' % (i)]
# #         del nparams['w%d' % (i)]
#         del nparams['tp%d' % (i)]
#     # new basis
    
    mod = radvel.RVModel(params, time_base=time_base)
    mod.params['dvdt'] = radvel.Parameter(value=0.0)
    mod.params['curv'] = radvel.Parameter(value=0.0)
    return mod


def make_like(mod,invels,planets,addextra=False):
    like = radvel.likelihood.RVLikelihood(mod, invels['time'], invels['mnvel'], invels['errvel'])
    gamma_guess = (invels['mnvel'].max() + invels['mnvel'].min())/2.
    
    like.params['gamma'] = radvel.Parameter(value=0.01)
    like.params['jit']= radvel.Parameter(value=3.)
    np = like.params.num_planets
    if addextra:
        np = np+1
    for n in range(1,np+1):
        like.params['w%d' % (n)].vary = True
        like.params['e%d'% (n)].vary = False
        like.params['per%d' % (n)].vary = False
        like.params['tp%d' % (n)].vary = False
#        like.params['secosw%d' % (n)].vary = False
#        like.params['sesinw%d'% (n)].vary = False
#        like.params['per%d' % (n)].vary = False
#        like.params['tc%d' % (n)].vary = True


        
    like.params['curv'].vary = False
    like.params['jit'].vary = True
    like.params['gamma'].vary = True        
    like.params['dvdt'].vary = False    
    
    return like

def init_posterior(like,planets,addextra=False):

    post = radvel.posterior.Posterior(like)
    post.priors += [radvel.prior.Gaussian( 'jit', 3, 3)]
    post.priors += [radvel.prior.Gaussian( 'gamma', float(like.params['gamma'].value), 30)]
    n = 1
    ntot = len(planets)
    if addextra:
        ntot = ntot + 1
        
#    for idx in range(0,ntot):
#        n = idx + 1
#        post.priors += [radvel.prior.Gaussian( 'k%d' % (n), np.log(5), 10)]
    post.priors += [radvel.prior.PositiveKPrior(ntot)]
    return post

def plot_results(post,chains,outdir,sname):

    radvel.plotting.corner_plot(post, chains,saveplot=os.path.join(outdir,sname+"_corner_cps_efix.pdf"))
    radvel.plotting.rv_multipanel_plot(post,saveplot=os.path.join(outdir,sname+"_mp_cps_efix.pdf"))

    return

def writefit(sname,outdir,periods,Ks,err_Ks,Ms,err_Ms,Kmls,rplanets,nvels):

    outname = "%s_efix.radvelfit" % (sname)
    outname = os.path.join(outdir,outname)
    outfp = open(outname,"w")
    for i in range(0,len(Ks)):
        ostr = "%f %f %f %f %f %f %d %f\n" % (periods[i],Ks[i],err_Ks[i],Ms[i],err_Ms[i],rplanets[i],nvels,Kmls[i])
        # period (days) K (m/s) err_K (m/s) planet mass (M_earth) error (M_earth) R (R_earth) #vs 
        outfp.write(ostr)

    outfp.close()

def mcmc_planets(post,outdir,sname,mstars,addextra=False):
    
    conf_base = sname
    nwalkers = 50
    nsteps = 10000
    ensembles = 8
    msg = "Running MCMC for {}, N_ensembles = {}, N_walkers = {}, N_steps = {} ...".format(
        conf_base, ensembles, nwalkers, nsteps)
    print msg

    chains = radvel.mcmc(
        post, nwalkers=nwalkers, nrun=nsteps,
        ensembles=ensembles
    )

    Ks = []
    err_Ks = []
    Ms = []
    err_Ms = []
    # Get quantiles and update posterior object
    post_summary=chains.quantile([0.159, 0.5, 0.841])
    c = 1
    if addextra:
        c = 2
    for n in range(1,post.params.num_planets+c):
        post_summary['Mpsini%d' % (n)] = radvel.utils.Msini((post_summary['k%d' % (n)]),post.params['per%d' % (n)].value,mstars[0],0.)
        if writefit:
#            Ks.append(post.params['k%d' % (n)].value)
            Ks.append(post_summary['k%d' % (n)][0.5])
            sigK = ( post_summary['k%d' % (n)][0.841] - post_summary['k%d' % (n)][0.159] ) / 2.
            sigK = sigK #* np.exp(post_summary['logk%d' % (n)][0.5])
            err_Ks.append(sigK)
#            Ms.append(radvel.utils.Msini((post.params['k%d' % (n)]).value,post.params['per%d' % (n)].value,mstars[0],0.))
            Ms.append(post_summary['Mpsini%d' % (n)][0.5])
            err_Ms.append(( post_summary['Mpsini%d' % (n)][0.841] - post_summary['Mpsini%d' % (n)][0.159]) / 2.)




    print "Saving output files..."
    saveto = os.path.join(outdir, sname+'_post_summary_cps_efix_wfix.csv')
    post_summary.to_csv(saveto, sep=',')
    
    return chains,Ks,err_Ks,Ms,err_Ms

if __name__ == "__main__":

    snames,gdfn,mfn,veldir,outdir,addextra = parse_options()

    gd = ascii.read(gdfn)
    TESSAPFdata=ascii.read(mfn,format='csv')
    # planet masses are in earth masses, Systemic likes Jupiters
    #mearth = 5.9722e24
    #mjup = 1.898e27
    mratio = 317.83
    TESSAPFdata['est_mass'] /= mratio
            
    for sname in snames:

        sfn = sname + ".sys"
        vfn = sname + ".vels"
        newvels = True

        invels = readin_velsfile(os.path.join(veldir,vfn))
#        ddates,dphases, dvels, derrs, di2sums = bin_phase_dates(invels["time"],invels["pphase"],invels['mnvel'],invels["errvel"],invels["i2val"])
#        bvfn = sname + "binned"
#        ascii.write([ddates,dvels,derrs,di2sums,dphases], os.path.join(veldir,bvfn+".vels"),format="no_header")
            
        if len(invels) > 4 :
            planets, = np.where((TESSAPFdata['star_names'] == sname) & (TESSAPFdata['detected'] == "TRUE"))

            mod = initialize_model(planets,TESSAPFdata,addextra=addextra)
            like = make_like(mod,invels,planets,addextra=addextra)
            post = init_posterior(like,planets,addextra=addextra)
            post = radvel.fitting.maxlike_fitting(post, verbose=False)
            print post
            for i in range(0,len(planets)):
                Kmls.append(post.params['k%d' % (i+1)].value)
            chains,Ks,err_Ks,Ms,err_Ms = mcmc_planets(post,outdir,sname,TESSAPFdata['mstar'][planets],addextra=addextra)
#            write_output(post,TESSAPFdata['phase'][planets],TESSAPFdata['vmag'][planets],TESSAPFdata['rplanet'][planets],TESSAPFdata['true_mass'][planets],TESSAPFdata['Index'][planets],gd['vel_offset'][gd['starname'] == sname],writefit=True,veldir=veldir,outdir=outdir)
            plot_results(post,chains,outdir,sname)
            writefit(sname,outdir,TESSAPFdata['period'][planets],Ks,err_Ks,Ms,err_Ms,Kmls,TESSAPFdata['rplanet'][planets],len(invels))
